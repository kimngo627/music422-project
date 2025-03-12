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
    N = 2*halfN
    shortHalfN = N_SHORT // 2
    shortN = N_SHORT
    # vectorizing the Dequantize function call
    #vDequantize = np.vectorize(Dequantize)

    if blockType == SHORT:
        num_short_blocks = (N // shortHalfN) - 1
        sfBands = ScaleFactorBands(AssignMDCTLinesFromFreqLimits(shortHalfN, codingParams.sampleRate, flimit=short_cbFreqLimits))

        data = np.zeros(N, dtype=np.float64)

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
            data[start:end] += SineWindow(IMDCT(mdctLine, shortHalfN, shortHalfN))

    else:
        rescaleLevel = 1.*(1<<overallScaleFactor)
        # reconstitute the first halfN MDCT lines of this channel from the stored data
        mdctLine = np.zeros(halfN,dtype=np.float64)
        iMant = 0
        for iBand in range(codingParams.sfBands.nBands):
            nLines =codingParams.sfBands.nLines[iBand]
            if bitAlloc[iBand]:
                mdctLine[iMant:(iMant+nLines)]=vDequantize(scaleFactor[iBand], mantissa[iMant:(iMant+nLines)],codingParams.nScaleBits, bitAlloc[iBand])
            iMant += nLines
        mdctLine /= rescaleLevel  # put overall gain back to original level


        # IMDCT and window the data for this channel
        if blockType == LONG:
            data = SineWindow( IMDCT(mdctLine, halfN, halfN) )  # takes in halfN MDCT coeffs
        elif blockType == START:
            data = TransitionWindow(IMDCT(mdctLine, halfN, shortHalfN), SineWindowFunc, N, shortN, start=True)
        elif blockType == STOP:
            data = TransitionWindow(IMDCT(mdctLine, halfN, shortHalfN), SineWindowFunc, N, shortN, start=False)

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
    
    timeSamples = data
    # window data for side chain FFT and also window and compute MDCT
    if blockType == LONG:
        mdctTimeSamples = SineWindow(data)
        mdctLines = MDCT(mdctTimeSamples, halfN, halfN)[:halfN]
        sfBands = codingParams.sfBands

    elif blockType == SHORT:
        num_short_blocks = (N // shortHalfN) - 1
        mdctLines = []
        
        for i in range(num_short_blocks):
            start = i * shortHalfN
            end = start + shortN
            short_data = data[start:end]
            short_windowed = SineWindow(short_data)
            short_mdctLines = MDCT(short_windowed, shortHalfN, shortHalfN)[:shortHalfN]
            mdctLines.append(short_mdctLines)

        sfBands = ScaleFactorBands(AssignMDCTLinesFromFreqLimits(shortHalfN, codingParams.sampleRate, flimit=short_cbFreqLimits))

    elif blockType == START:
        mdctTimeSamples = TransitionWindow(data, SineWindowFunc, N, shortN, start=True)
        mdctLines = MDCT(mdctTimeSamples, halfN, shortHalfN)[:halfN]
        sfBands = codingParams.sfBands

    elif blockType == STOP:
        mdctTimeSamples = TransitionWindow(data, SineWindowFunc, N, shortN, start=False)
        mdctLines = MDCT(mdctTimeSamples, shortHalfN, halfN)[:halfN]
        sfBands = codingParams.sfBands

    # compute target mantissa bit budget for this block of halfN MDCT mantissas
    # targetBitsPerSample also changes because we took out some bits for the header
    bitBudget = codingParams.targetBitsPerSample * halfN  # this is overall target bit rate    # TODO may need to change targetBitsPerSample
    bitBudget -=  nScaleBits*(sfBands.nBands +1)  # less scale factor bits (including overall scale factor)
    bitBudget -= codingParams.nMantSizeBits*sfBands.nBands  # less mantissa bit allocation bits

    if blockType == SHORT: 
        scaleFactor = []
        bitAlloc = []
        mantissa = []
        overallScale = []

        for i in range(num_short_blocks):
            start = i * shortHalfN
            end = start + shortN
            short_mdctLines = mdctLines[i]
            
            # compute overall scale factor for this block and boost mdctLines using it
            maxLine = np.max( np.abs(short_mdctLines) )
            short_overallScale = ScaleFactor(maxLine,nScaleBits)  #leading zeroes don't depend on nMantBits
            short_mdctLines *= (1<<short_overallScale)

            # compute the mantissa bit allocations
            # compute SMRs in side chain FFT
            SMRs = CalcSMRs(timeSamples[start:end], short_mdctLines, overallScale, codingParams.sampleRate, sfBands)
            # perform bit allocation using SMR results
            short_bitAlloc = BitAlloc(bitBudget//num_short_blocks, maxMantBits, sfBands.nBands, sfBands.nLines, SMRs)

            # given the bit allocations, quantize the mdct lines in each band
            short_scaleFactor = np.empty(sfBands.nBands,dtype=np.int32)
            nMant=shortHalfN
            for iBand in range(sfBands.nBands):
                if not short_bitAlloc[iBand]: nMant-= sfBands.nLines[iBand]  # account for mantissas not being transmitted
            short_mantissa=np.empty(nMant,dtype=np.int32)
            iMant=0
            for iBand in range(sfBands.nBands):
                lowLine = sfBands.lowerLine[iBand]
                highLine = sfBands.upperLine[iBand] + 1  # extra value is because slices don't include last value
                nLines= sfBands.nLines[iBand]
                scaleLine = np.max(np.abs( short_mdctLines[lowLine:highLine] ) )
                short_scaleFactor[iBand] = ScaleFactor(scaleLine, nScaleBits, short_bitAlloc[iBand])
                if short_bitAlloc[iBand]:
                    short_mantissa[iMant:iMant+nLines] = vMantissa(short_mdctLines[lowLine:highLine],short_scaleFactor[iBand], nScaleBits, short_bitAlloc[iBand])
                    iMant += nLines
            # end of loop over scale factor bands

            scaleFactor.append(short_scaleFactor)
            bitAlloc.append(short_bitAlloc)
            mantissa.append(short_mantissa)
            overallScale.append(short_overallScale)
    
    else:

        # compute overall scale factor for this block and boost mdctLines using it
        maxLine = np.max( np.abs(mdctLines) )
        overallScale = ScaleFactor(maxLine,nScaleBits)  #leading zeroes don't depend on nMantBits
        mdctLines *= (1<<overallScale)

        # compute the mantissa bit allocations
        # compute SMRs in side chain FFT
        SMRs = CalcSMRs(timeSamples, mdctLines, overallScale, codingParams.sampleRate, sfBands)
        # perform bit allocation using SMR results
        bitAlloc = BitAlloc(bitBudget, maxMantBits, sfBands.nBands, sfBands.nLines, SMRs)

        # given the bit allocations, quantize the mdct lines in each band
        scaleFactor = np.empty(sfBands.nBands,dtype=np.int32)
        nMant=halfN
        for iBand in range(sfBands.nBands):
            if not bitAlloc[iBand]: nMant-= sfBands.nLines[iBand]  # account for mantissas not being transmitted
        mantissa=np.empty(nMant,dtype=np.int32)
        iMant=0
        for iBand in range(sfBands.nBands):
            lowLine = sfBands.lowerLine[iBand]
            highLine = sfBands.upperLine[iBand] + 1  # extra value is because slices don't include last value
            nLines= sfBands.nLines[iBand]
            scaleLine = np.max(np.abs( mdctLines[lowLine:highLine] ) )
            scaleFactor[iBand] = ScaleFactor(scaleLine, nScaleBits, bitAlloc[iBand])
            if bitAlloc[iBand]:
                mantissa[iMant:iMant+nLines] = vMantissa(mdctLines[lowLine:highLine],scaleFactor[iBand], nScaleBits, bitAlloc[iBand])
                iMant += nLines
        # end of loop over scale factor bands

    # return results
    return (scaleFactor, bitAlloc, mantissa, overallScale)



