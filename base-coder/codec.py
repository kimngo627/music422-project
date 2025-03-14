"""
codec.py -- The actual encode/decode functions for the perceptual audio codec

-----------------------------------------------------------------------
© 2019-2025 Marina Bosi & Richard E. Goldberg -- All rights reserved
-----------------------------------------------------------------------
"""

import numpy as np  # used for arrays

# used by Encode and Decode
from window import SineWindow  # current window used for MDCT -- implement KB-derived?
from mdct import MDCT,IMDCT  # fast MDCT implementation (uses numpy FFT)
from quantize import *  # using vectorized versions (to use normal versions, uncomment lines 18,67 below defining vMantissa and vDequantize)

# used only by Encode
from psychoac import *  # calculates SMRs for each scale factor band
from bitalloc import BitAlloc  #allocates bits to scale factor bands given SMRs
from transient_detector import BlockType
from blockswitching import *

# binary used for keeping track of what kind of block should be used
LONG  = bin(0)
START = bin(1)
SHORT = bin(2)
STOP  = bin(3)

N_SHORT = 128

def Decode(scaleFactor,bitAlloc,mantissa,overallScaleFactor,codingParams):
    """Reconstitutes a single-channel block of encoded data into a block of
    signed-fraction data based on the parameters in a PACFile object"""
    blockType = codingParams.currBlockType
    halfN = codingParams.nMDCTLines
    N = halfN*2
    shortHalfN = N_SHORT // 2
    shortN = N_SHORT
    # vectorizing the Dequantize function call
    #vDequantize = np.vectorize(Dequantize)
    data = np.zeros(N, dtype=np.float64)

    if blockType == SHORT:
        num_short_blocks = (halfN // N_SHORT) * 2 - 1
        if hasattr(codingParams, 'shortSfBands'):
            sfBands = codingParams.shortSfBands
        else:
            # Create them if they don't exist
            sfBands = ScaleFactorBands(
                AssignMDCTLinesFromFreqLimits(shortHalfN, codingParams.sampleRate, 
                                             flimit=short_cbFreqLimits)
            )
        fullShortData = np.zeros(N, dtype=np.float64)
        for i in range(num_short_blocks):
            rescaleLevel = 1.*(1<<overallScaleFactor[i])
            mdctLine = np.zeros(shortHalfN, dtype=np.float64)
            iMant = 0

            for iBand in range(sfBands.nBands):
                nLines = sfBands.nLines[iBand]
                if bitAlloc[i][iBand]:
                    mdctLine[iMant:(iMant+nLines)] = vDequantize(scaleFactor[i][iBand], mantissa[i][iMant:(iMant+nLines)],codingParams.nScaleBits, bitAlloc[i][iBand])
                iMant += nLines
            mdctLine /= rescaleLevel

            start = i * shortHalfN
            end = start + shortN
            if end > N: end = N

            short_data = SineWindow(IMDCT(mdctLine, shortHalfN, shortHalfN))
            fullShortData[start:end] += short_data[:shortN]
        return fullShortData

    else:
        rescaleLevel = 1.*(1<<overallScaleFactor)
        # reconstitute the first halfN MDCT lines of this channel from the stored data
        if blockType == LONG:
            mdctLine = np.zeros(halfN,dtype=np.float64)
            sfBands = codingParams.sfBands
        elif blockType == START:
            mdctLine = np.zeros((halfN + shortHalfN)//2,dtype=np.float64)
            sfBands = codingParams.transitionSfBands
        elif blockType == STOP:
            mdctLine = np.zeros((halfN + shortHalfN)//2,dtype=np.float64)
            sfBands = codingParams.transitionSfBands
        iMant = 0
        for iBand in range(sfBands.nBands):
            nLines =sfBands.nLines[iBand]
            if bitAlloc[iBand]:
                mdctLine[iMant:(iMant+nLines)]=vDequantize(scaleFactor[iBand], mantissa[iMant:(iMant+nLines)],codingParams.nScaleBits, bitAlloc[iBand])
            iMant += nLines
        mdctLine /= rescaleLevel  # put overall gain back to original level


        # IMDCT and window the data for this channel
        if blockType == LONG:
            data = SineWindow( IMDCT(mdctLine, halfN, halfN) )  # takes in halfN MDCT coeffs
        elif blockType == START:
            # need to use halfN and shortHalfN
            data = TransitionWindow(IMDCT(mdctLine, halfN, shortHalfN), SineWindowFunc, N, shortN, start=True)
        elif blockType == STOP:
            # need to use halfN and shortHalfN
            data = TransitionWindow(IMDCT(mdctLine, shortHalfN, halfN), SineWindowFunc, N, shortN, start=False)
        else:
            # Fallback for unexpected block types
            print(f"Warning: Unknown block type {blockType} in Decode, defaulting to LONG")
            data = SineWindow(IMDCT(mdctLine, halfN, halfN))

    # end loop over channels, return reconstituted time samples (pre-overlap-and-add)
        return data


def Encode(data,codingParams):
    """Encodes a multi-channel block of signed-fraction data based on the parameters in a PACFile object"""
    scaleFactor = []
    bitAlloc = []
    mantissa = []
    overallScaleFactor = []

    # loop over channels and separately encode each one
    for iCh in range(codingParams.nChannels):
        (s,b,m,o) = EncodeSingleChannel(data[iCh],codingParams)
        scaleFactor.append(s)
        bitAlloc.append(b)
        mantissa.append(m)
        overallScaleFactor.append(o)
    # return results bundled over channels
    return (scaleFactor,bitAlloc,mantissa,overallScaleFactor)


def EncodeSingleChannel(data,codingParams):
    """Encodes a single-channel block of signed-fraction data based on the parameters in a PACFile object"""

    # prepare various constants 
    halfN = codingParams.nMDCTLines
    N = 2*halfN
    shortHalfN = N_SHORT // 2
    shortN = N_SHORT
    nScaleBits = codingParams.nScaleBits
    maxMantBits = (1<<codingParams.nMantSizeBits)  # 1 isn't an allowed bit allocation so n size bits counts up to 2^n
    if maxMantBits>16: maxMantBits = 16  # to make sure we don't ever overflow mantissa holders
    blockType = codingParams.currBlockType
    # vectorizing the Mantissa function call
    #vMantissa = np.vectorize(Mantissa)
    # TODO: define sfBands for transition blocks
    #mdctLines = MDCT(mdctTimeSamples, halfN, halfN)[:halfN]
    
    timeSamples = data
    # window data for side chain FFT and also window and compute MDCT
    if blockType == LONG:
        mdctTimeSamples = SineWindow(data)
        mdctLines = MDCT(mdctTimeSamples, halfN, halfN)[:halfN]
        sfBands = codingParams.sfBands

    elif blockType == SHORT:
        num_short_blocks = (halfN // N_SHORT) * 2 - 1
        mdctLines = []

        if hasattr(codingParams, 'shortSfBands'):
            sfBands = codingParams.shortSfBands
        else:
            sfBands = ScaleFactorBands(
                AssignMDCTLinesFromFreqLimits(shortHalfN, codingParams.sampleRate, 
                                             flimit=short_cbFreqLimits)
            )
        
        for i in range(num_short_blocks):
            start = i * shortHalfN
            end = start + shortN
            short_data = data[start:end]
            short_windowed = SineWindow(short_data)
            short_mdctLines = MDCT(short_windowed, shortHalfN, shortHalfN)[:shortHalfN]
            mdctLines.append(short_mdctLines)

    elif blockType == START:
        mdctTimeSamples = TransitionWindow(data, SineWindowFunc, N, shortN, start=True)
        mdctLines = MDCT(mdctTimeSamples, halfN, shortHalfN)[:halfN]
        sfBands = codingParams.transitionSfBands

    elif blockType == STOP:
        mdctTimeSamples = TransitionWindow(data, SineWindowFunc, N, shortN, start=False)
        mdctLines = MDCT(mdctTimeSamples, shortHalfN, halfN)[:halfN]
        sfBands = codingParams.transitionSfBands

    else:
        # Fallback for undefined block types (this should never happen)
        print(f"Warning: Unknown block type {blockType}, defaulting to LONG")
        mdctTimeSamples = SineWindow(data)
        mdctLines = MDCT(mdctTimeSamples, halfN, halfN)[:halfN]

    # compute target mantissa bit budget for this block of halfN MDCT mantissas
    # targetBitsPerSample also changes because we took out some bits for the header
    bitBudget = codingParams.targetBitsPerSample * halfN  # this is overall target bit rate    

    if blockType == SHORT:
        num_short_blocks = (halfN // N_SHORT) * 2 - 1
        sfBands = codingParams.shortSfBands
        
        bitBudget -= nScaleBits * num_short_blocks 
        bitBudget -= (nScaleBits + codingParams.nMantSizeBits) * sfBands.nBands * num_short_blocks
        
        scaleFactor = []
        bitAlloc = []
        mantissa = []
        overallScale = []

        for i in range(num_short_blocks):
            start = i * shortHalfN
            end = start + shortN
            short_mdctLines = mdctLines[i]
            
            maxLine = np.max(np.abs(short_mdctLines))
            short_overallScale = ScaleFactor(maxLine, nScaleBits)
            short_mdctLines *= (1 << short_overallScale)

            analysis_data = timeSamples[start:end]
            if len(analysis_data) < shortN:
                analysis_data = np.pad(analysis_data, (0, shortN - len(analysis_data)))

            SMRs = CalcSMRs(
                analysis_data, 
                short_mdctLines, 
                short_overallScale,  
                codingParams.sampleRate, 
                sfBands
            )
            
            short_bitBudget = bitBudget // num_short_blocks
            short_bitAlloc = BitAlloc(short_bitBudget, maxMantBits, sfBands.nBands, 
                                      sfBands.nLines, SMRs)

            short_scaleFactor = np.empty(sfBands.nBands, dtype=np.int32)
            
            nMant = shortHalfN
            for iBand in range(sfBands.nBands):
                if not short_bitAlloc[iBand]: 
                    nMant -= sfBands.nLines[iBand]
                    
            short_mantissa = np.empty(nMant, dtype=np.int32)
            iMant = 0
            
            for iBand in range(sfBands.nBands):
                lowLine = sfBands.lowerLine[iBand]
                highLine = sfBands.upperLine[iBand] + 1
                nLines = sfBands.nLines[iBand]
                
                scaleLine = np.max(np.abs(short_mdctLines[lowLine:highLine]))
                short_scaleFactor[iBand] = ScaleFactor(scaleLine, nScaleBits, short_bitAlloc[iBand])
                
                if short_bitAlloc[iBand]:
                    short_mantissa[iMant:iMant+nLines] = vMantissa(
                        short_mdctLines[lowLine:highLine],
                        short_scaleFactor[iBand], 
                        nScaleBits, 
                        short_bitAlloc[iBand]
                    )
                    iMant += nLines

            scaleFactor.append(short_scaleFactor)
            bitAlloc.append(short_bitAlloc)
            mantissa.append(short_mantissa)
            overallScale.append(short_overallScale)
    
    else:
        # For LONG, START, STOP blocks
        if codingParams.currBlockType == LONG:
            sfBands = codingParams.sfBands
        else:
            sfBands = codingParams.transitionSfBands
        bitBudget -= nScaleBits  
        bitBudget -= (nScaleBits + codingParams.nMantSizeBits) * sfBands.nBands

        maxLine = np.max(np.abs(mdctLines))
        overallScale = ScaleFactor(maxLine, nScaleBits)
        mdctLines *= (1 << overallScale)

        required_length = sfBands.upperLine[-1] + 1
        if len(mdctLines) < required_length:
            # Pad mdctLines if it's too short
            mdctLines = np.pad(mdctLines, (0, required_length - len(mdctLines)))


        SMRs = CalcSMRs(timeSamples, mdctLines, overallScale, codingParams.sampleRate, sfBands)
        bitAlloc = BitAlloc(bitBudget, maxMantBits, sfBands.nBands, sfBands.nLines, SMRs)

        scaleFactor = np.empty(sfBands.nBands, dtype=np.int32)

        nMant = halfN
        for iBand in range(sfBands.nBands):
            if not bitAlloc[iBand]: 
                nMant -= sfBands.nLines[iBand]
                
        mantissa = np.empty(nMant, dtype=np.int32)
        iMant = 0
        
        for iBand in range(sfBands.nBands):
            lowLine = sfBands.lowerLine[iBand]
            highLine = sfBands.upperLine[iBand] + 1
            nLines = sfBands.nLines[iBand]
            
            scaleLine = np.max(np.abs(mdctLines[lowLine:highLine]))
            scaleFactor[iBand] = ScaleFactor(scaleLine, nScaleBits, bitAlloc[iBand])
            
            if bitAlloc[iBand]:
                mantissa[iMant:iMant+nLines] = vMantissa(
                    mdctLines[lowLine:highLine],
                    scaleFactor[iBand], 
                    nScaleBits, 
                    bitAlloc[iBand]
                )
                iMant += nLines

    # return results
    return (scaleFactor, bitAlloc, mantissa, overallScale)



