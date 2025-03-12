"""
pacfile.py -- Defines a PACFile class to handle reading and writing audio
data to an audio file holding data compressed using an MDCT-based perceptual audio
coding algorithm.  The MDCT lines of each audio channel are grouped into bands,
each sharing a single scaleFactor and bit allocation that are used to block-
floating point quantize those lines.  This class is a subclass of AudioFile.

-----------------------------------------------------------------------
© 2019-2025 Marina Bosi & Richard E. Goldberg -- All rights reserved
-----------------------------------------------------------------------

See the documentation of the AudioFile class for general use of the AudioFile
class.

Notes on reading and decoding PAC files:

    The OpenFileForReading() function returns a CodedParams object containing:

        nChannels = the number of audio channels
        sampleRate = the sample rate of the audio samples
        numSamples = the total number of samples in the file for each channel
        nMDCTLines = half the MDCT block size (block switching not supported)
        nSamplesPerBlock = MDCTLines (but a name that PCM files look for)
        nScaleBits = the number of bits storing scale factors
        nMantSizeBits = the number of bits storing mantissa bit allocations
        sfBands = a ScaleFactorBands object
        overlapAndAdd = decoded data from the prior block (initially all zeros)

    The returned ScaleFactorBands object, sfBands, contains an allocation of
    the MDCT lines into groups that share a single scale factor and mantissa bit
    allocation.  sfBands has the following attributes available:

        nBands = the total number of scale factor bands
        nLines[iBand] = the number of MDCT lines in scale factor band iBand
        lowerLine[iBand] = the first MDCT line in scale factor band iBand
        upperLine[iBand] = the last MDCT line in scale factor band iBand


Notes on encoding and writing PAC files:

    When writing to a PACFile the CodingParams object passed to OpenForWriting()
    should have the following attributes set:

        nChannels = the number of audio channels
        sampleRate = the sample rate of the audio samples
        numSamples = the total number of samples in the file for each channel
        nMDCTLines = half the MDCT block size (format does not support block switching)
        nSamplesPerBlock = MDCTLines (but a name that PCM files look for)
        nScaleBits = the number of bits storing scale factors
        nMantSizeBits = the number of bits storing mantissa bit allocations
        targetBitsPerSample = the target encoding bit rate in units of bits per sample

    The first three attributes (nChannels, sampleRate, and numSamples) are
    typically added by the original data source (e.g. a PCMFile object) but
    numSamples may need to be extended to account for the MDCT coding delay of
    nMDCTLines and any zero-padding done in the final data block

    OpenForWriting() will add the following attributes to be used during the encoding
    process carried out in WriteDataBlock():

        sfBands = a ScaleFactorBands object
        priorBlock = the prior block of audio data (initially all zeros)

    The passed ScaleFactorBands object, sfBands, contains an allocation of
    the MDCT lines into groups that share a single scale factor and mantissa bit
    allocation.  sfBands has the following attributes available:

        nBands = the total number of scale factor bands
        nLines[iBand] = the number of MDCT lines in scale factor band iBand
        lowerLine[iBand] = the first MDCT line in scale factor band iBand
        upperLine[iBand] = the last MDCT line in scale factor band iBand

Description of the PAC File Format:

    Header:

        tag                 4 byte file tag equal to "PAC "
        sampleRate          little-endian unsigned long ("<L" format in struct)
        nChannels           little-endian unsigned short("<H" format in struct)
        numSamples          little-endian unsigned long ("<L" format in struct)
        nMDCTLines          little-endian unsigned long ("<L" format in struct)
        nScaleBits          little-endian unsigned short("<H" format in struct)
        nMantSizeBits       little-endian unsigned short("<H" format in struct)
        nSFBands            little-endian unsigned long ("<L" format in struct)
        for iBand in range(nSFBands):
            nLines[iBand]   little-endian unsigned short("<H" format in struct)

    Each Data Block:  (reads data blocks until end of file hit)

        for iCh in range(nChannels):
            nBytes          little-endian unsigned long ("<L" format in struct)
            as bits packed into an array of nBytes bytes:
                overallScale[iCh]                       nScaleBits bits
                for iBand in range(nSFBands):
                    scaleFactor[iCh][iBand]             nScaleBits bits
                    bitAlloc[iCh][iBand]                nMantSizeBits bits
                    if bitAlloc[iCh][iBand]:
                        for m in nLines[iBand]:
                            mantissa[iCh][iBand][m]     bitAlloc[iCh][iBand]+1 bits
                <extra custom data bits as long as space is included in nBytes>

"""
'''
Changes from block switching 
ReadFileHeader(): 
    myParams.numSamples = numSamples + nMDCTLines   # because we prepend N zeros to the beginning of audio file
WriteDataBlock(): additional coding parameters 
    # binary used for keeping track of what kind of block should be used
    # 0b00 = regular long, 0b01 = start, 0b10 = short,  0b11 = stop
    codingParams.currBlockType = bin(0)
    # we are using 8 short blocks in substitution of a long block, keep track which one
    codingParams.shortBlockInd = bin(7)  # indexed 0-7
    # codingParams.shortBlockInd used to indicate if we should use start or stop block

'''

from audiofile import * # base class
from bitpack import *  # class for packing data into an array of bytes where each item's number of bits is specified
import codec    # module where the actual PAC coding functions reside(this module only specifies the PAC file format)
from psychoac import ScaleFactorBands, AssignMDCTLinesFromFreqLimits  # defines the grouping of MDCT lines into scale factor bands

import numpy as np  # to allow conversion of data blocks to numpy's array object
MAX16BITS = 32767

from transient_detector import *
SR = 44100

# binary used for keeping track of what kind of block should be used
LONG  = bin(0)
START = bin(1)
SHORT = bin(2)
STOP  = bin(3)

N = 1024

class PACFile(AudioFile):
    """
    Handlers for a perceptually coded audio file I am encoding/decoding
    """

    # a file tag to recognize PAC coded files
    tag=b'PAC '

    def ReadFileHeader(self):
        """
        Reads the PAC file header from a just-opened PAC file and uses it to set
        object attributes.  File pointer ends at start of data portion.
        """
        # check file header tag to make sure it is the right kind of file
        tag=self.fp.read(4)
        if tag!=self.tag: raise RuntimeError("Tried to read a non-PAC file into a PACFile object")
        # use struct.unpack() to load up all the header data
        (sampleRate, nChannels, numSamples, nMDCTLines, nScaleBits, nMantSizeBits) \
                 = unpack('<LHLLHH',self.fp.read(calcsize('<LHLLHH')))
        nBands = unpack('<L',self.fp.read(calcsize('<L')))[0]
        nLines=  unpack('<'+str(nBands)+'H',self.fp.read(calcsize('<'+str(nBands)+'H')))
        sfBands=ScaleFactorBands(nLines)
        # load up a CodingParams object with the header data
        myParams=CodingParams()
        myParams.sampleRate = sampleRate
        myParams.nChannels = nChannels
        myParams.numSamples = numSamples # + nMDCTLines   # because we prepend N zeros to the beginning of audio file
        myParams.nMDCTLines = myParams.nSamplesPerBlock = nMDCTLines
        myParams.nScaleBits = nScaleBits
        myParams.nMantSizeBits = nMantSizeBits
        # add in scale factor band information
        myParams.sfBands =sfBands
        # start w/o all zeroes as data from prior block to overlap-and-add for output
        overlapAndAdd = []
        for iCh in range(nChannels): overlapAndAdd.append( np.zeros(nMDCTLines, dtype=np.float64) )
        myParams.overlapAndAdd=overlapAndAdd
        return myParams


    def ReadDataBlock(self, codingParams):
        """
        Reads a block of coded data from a PACFile object that has already
        executed OpenForReading() and returns those samples as reconstituted
        signed-fraction data
        """
        # loop over channels (whose coded data are stored separately) and read in each data block
        data=[]
        for iCh in range(codingParams.nChannels):
            data.append(np.array([],dtype=np.float64))  # add location for this channel's data
            # read in string containing the number of bytes of data for this channel (but check if at end of file!)
            s=self.fp.read(calcsize("<L"))  # will be empty if at end of file
            if not s:
                # hit last block, see if final overlap and add needs returning, else return nothing
                if codingParams.overlapAndAdd:
                    overlapAndAdd=codingParams.overlapAndAdd
                    codingParams.overlapAndAdd=0  # setting it to zero so next pass will just return
                    return overlapAndAdd
                else:
                    return
            # not at end of file, get nBytes from the string we just read
            nBytes = unpack("<L",s)[0] # read it as a little-endian unsigned long
            # read the nBytes of data into a PackedBits object to unpack
            pb = PackedBits()
            pb.SetPackedData( self.fp.read(nBytes) ) # PackedBits function SetPackedData() converts strings to internally-held array of bytes
            if pb.nBytes < nBytes:  raise "Only read a partial block of coded PACFile data"

            # Read block type (2 bits)
            blockType_bits = pb.ReadBits(2)
            if blockType_bits == 0:
                codingParams.currBlockType = LONG    # 00 = LONG
            elif blockType_bits == 1:
                codingParams.currBlockType = SHORT   # 01 = SHORT
            elif blockType_bits == 2:
                codingParams.currBlockType = START   # 10 = START
            elif blockType_bits == 3:
                codingParams.currBlockType = STOP    # 11 = STOP
            # print("block type: ", codingParams.currBlockType)
            # extract the data from the PackedBits object
            if codingParams.currBlockType == SHORT:
                # TODO: Modify the following code to loop through and read short blocks one by one
                num_short_blocks = (len(data) // (N_SHORT//2)) - 1
                for i in range(num_short_blocks):
                    overallScaleFactor = pb.ReadBits(codingParams.nScaleBits)  # overall scale factor
                    scaleFactor=[]
                    bitAlloc=[]
                    mantissa=np.zeros(codingParams.nMDCTLines,np.int32)  # start w/ all mantissas zero
                    for iBand in range(codingParams.sfBands.nBands): # loop over each scale factor band to pack its data
                        ba = pb.ReadBits(codingParams.nMantSizeBits)
                        if ba: ba+=1  # no bit allocation of 1 so ba of 2 and up stored as one less
                        bitAlloc.append(ba)  # bit allocation for this band
                        scaleFactor.append(pb.ReadBits(codingParams.nScaleBits))  # scale factor for this band
                        if bitAlloc[iBand]:
                            # if bits allocated, extract those mantissas and put in correct location in matnissa array
                            m=np.empty(codingParams.sfBands.nLines[iBand],np.int32)
                            for j in range(codingParams.sfBands.nLines[iBand]):
                                m[j]=pb.ReadBits(bitAlloc[iBand])     # mantissas for this band (if bit allocation non-zero) and bit alloc <>1 so encoded as 1 lower than actual allocation
                            mantissa[codingParams.sfBands.lowerLine[iBand]:(codingParams.sfBands.upperLine[iBand]+1)] = m
                    # done unpacking data (end loop over scale factor bands)

                    # CUSTOM DATA:
                    # < now can unpack any custom data passed in the nBytes of data > 
                    # block type
                
                
                    # (DECODE HERE) decode the unpacked data for this channel, overlap-and-add first half, and append it to the data array (saving other half for next overlap-and-add)
                    decodedData = self.Decode(scaleFactor,bitAlloc,mantissa, overallScaleFactor,codingParams)
                    data[iCh] = np.concatenate( (data[iCh],np.add(codingParams.overlapAndAdd[iCh],decodedData[:codingParams.nMDCTLines]) ) )  # data[iCh] is overlap-and-added data
                    codingParams.overlapAndAdd[iCh] = decodedData[codingParams.nMDCTLines:]  # save other half for next pass

            else:
                overallScaleFactor = pb.ReadBits(codingParams.nScaleBits)  # overall scale factor
                scaleFactor=[]
                bitAlloc=[]
                mantissa=np.zeros(codingParams.nMDCTLines,np.int32)  # start w/ all mantissas zero
                for iBand in range(codingParams.sfBands.nBands): # loop over each scale factor band to pack its data
                    ba = pb.ReadBits(codingParams.nMantSizeBits)
                    if ba: ba+=1  # no bit allocation of 1 so ba of 2 and up stored as one less
                    bitAlloc.append(ba)  # bit allocation for this band
                    scaleFactor.append(pb.ReadBits(codingParams.nScaleBits))  # scale factor for this band
                    if bitAlloc[iBand]:
                        # if bits allocated, extract those mantissas and put in correct location in matnissa array
                        m=np.empty(codingParams.sfBands.nLines[iBand],np.int32)
                        for j in range(codingParams.sfBands.nLines[iBand]):
                            m[j]=pb.ReadBits(bitAlloc[iBand])     # mantissas for this band (if bit allocation non-zero) and bit alloc <>1 so encoded as 1 lower than actual allocation
                        mantissa[codingParams.sfBands.lowerLine[iBand]:(codingParams.sfBands.upperLine[iBand]+1)] = m
                # done unpacking data (end loop over scale factor bands)

                # CUSTOM DATA:
                # < now can unpack any custom data passed in the nBytes of data > 
                # block type
                
                # (DECODE HERE) decode the unpacked data for this channel, overlap-and-add first half, and append it to the data array (saving other half for next overlap-and-add)
                decodedData = self.Decode(scaleFactor,bitAlloc,mantissa, overallScaleFactor,codingParams)
                data[iCh] = np.concatenate( (data[iCh],np.add(codingParams.overlapAndAdd[iCh],decodedData[:codingParams.nMDCTLines]) ) )  # data[iCh] is overlap-and-added data
                codingParams.overlapAndAdd[iCh] = decodedData[codingParams.nMDCTLines:]  # save other half for next pass

        # end loop over channels, return signed-fraction samples for this block
        return data


    def WriteFileHeader(self,codingParams):
        """
        Writes the PAC file header for a just-opened PAC file and uses codingParams
        attributes for the header data.  File pointer ends at start of data portion.
        """
        # write a header tag
        self.fp.write(self.tag)
        # make sure that the number of samples in the file is a multiple of the
        # number of MDCT half-blocksize, otherwise zero pad as needed
        if not codingParams.numSamples%codingParams.nMDCTLines:
            codingParams.numSamples += (codingParams.nMDCTLines
                        - codingParams.numSamples%codingParams.nMDCTLines) # zero padding for partial final PCM block
        # also add in the delay block for the second pass w/ the last half-block
        codingParams.numSamples+= codingParams.nMDCTLines  # due to the delay in processing the first samples on both sides of the MDCT block
        # write the coded file attributes
        self.fp.write(pack('<LHLLHH',
            codingParams.sampleRate, codingParams.nChannels,
            codingParams.numSamples, codingParams.nMDCTLines,
            codingParams.nScaleBits, codingParams.nMantSizeBits  ))
        # create a ScaleFactorBand object to be used by the encoding process and write its info to header
        sfBands=ScaleFactorBands( AssignMDCTLinesFromFreqLimits(codingParams.nMDCTLines,    #TODO change number of mdctlines?
                                                                codingParams.sampleRate)
                                )
        # codingParams.secretMessage = "writing file header secret message"
        codingParams.sfBands=sfBands
        self.fp.write(pack('<L',sfBands.nBands))
        self.fp.write(pack('<'+str(sfBands.nBands)+'H',*(sfBands.nLines.tolist()) ))
        # start w/o all zeroes as prior block of unencoded data for other half of MDCT block
        priorBlock = []
        for iCh in range(codingParams.nChannels):
            priorBlock.append(np.zeros(codingParams.nMDCTLines,dtype=np.float64) )
        codingParams.priorBlock = priorBlock
        codingParams.currBlockType = LONG
        return

    ''''
        1. prepend N zeros infront of first block, enter while loop. -- already included in starter code 
        2. readdatablock
        3. writedatablock:
            1. we will detect for transience in first block (current block). If find transience, 
            block type of zeros is start block. current block type is short. 
            2. do mdct on the zeros + first half of first block.
        - next iteration of while loop. repeat #1,2,3
            1. we will detect for transience in second block.
            2. do mdct on N/2 zeros + first half of first block.
        - next iteration of while loop. repeat #1,2,3, 3.1, 3.2
    '''
    def WriteDataBlock(self,data, codingParams):
        """
        Writes a block of signed-fraction data to a PACFile object that has
        already executed OpenForWriting()"""

        # combine this block of multi-channel data w/ the prior block's to prepare for MDCTs twice as long
        fullBlockData=[]
        for iCh in range(codingParams.nChannels):
            fullBlockData.append( np.concatenate( ( codingParams.priorBlock[iCh], data[iCh]) ) )
        codingParams.priorBlock = data  # current pass's data is next pass's prior block data
        
        # print("inside write data block ", codingParams.__dict__)
        # we want to first start off detecting for transience in next block
        next_block = inFile.ReadDataBlock(codingParams)
        # print("call to read data block inside write data block ", codingParams.__dict__)
        if next_block:
            # print(len(next_block[0]))
            is_transient = detector.detect(next_block[0], SR)

            if hasattr(codingParams, 'blockState'):
                next_block_type = codingParams.blockState.update(is_transient)
            else:
                if codingParams.currBlockType == LONG and is_transient:
                    next_block_type = START
                elif codingParams.currBlockType == START:
                    next_block_type = SHORT
                elif codingParams.currBlockType == SHORT and is_transient:
                    next_block_type = SHORT
                elif codingParams.currBlockType == SHORT and not is_transient:
                    next_block_type = STOP
                elif codingParams.currBlockType == STOP:
                    next_block_type = LONG
                else:
                    next_block_type = LONG

            codingParams.nextBlockType = next_block_type
        else:
            codingParams.nextBlockType = LONG

        # (ENCODE HERE) Encode the full block of multi=channel data
        (scaleFactor,bitAlloc,mantissa, overallScaleFactor) = self.Encode(fullBlockData,codingParams)  # returns a tuple with all the block-specific info not in the file header

        # for each channel, write the data to the output file
        for iCh in range(codingParams.nChannels):

            # determine the size of this channel's data block and write it to the output file
            nBytes =codingParams.nScaleBits  # bits for overall scale factor
            # Add 2 bits for block type information
            nBytes += 2

            if codingParams.currBlockType == SHORT:
                num_short_blocks = 15
                for i in range(num_short_blocks):
                    for iBand in range(codingParams.sfBands.nBands): # loop over each scale factor band to get its bits
                        nBytes += codingParams.nMantSizeBits+codingParams.nScaleBits    # mantissa bit allocation and scale factor for that sf band
                        if bitAlloc[iCh][i][iBand]:
                            # if non-zero bit allocation for this band, add in bits for scale factor and each mantissa (0 bits means zero)
                            nBytes += bitAlloc[iCh][i][iBand]*codingParams.sfBands.nLines[iBand]  # no bit alloc = 1 so actuall alloc is one higher
            else:
                for iBand in range(codingParams.sfBands.nBands): # loop over each scale factor band to get its bits
                    nBytes += codingParams.nMantSizeBits+codingParams.nScaleBits    # mantissa bit allocation and scale factor for that sf band
                    if bitAlloc[iCh][iBand]:
                        # if non-zero bit allocation for this band, add in bits for scale factor and each mantissa (0 bits means zero)
                        nBytes += bitAlloc[iCh][iBand]*codingParams.sfBands.nLines[iBand]  # no bit alloc = 1 so actuall alloc is one higher
            # end computing bits needed for this channel's data

            # CUSTOM DATA:
            # < now can add space for custom data, if desired> # TODO for specific instructions on block? 

            # now convert the bits to bytes (w/ extra one if spillover beyond byte boundary)
            if nBytes%BYTESIZE==0:  nBytes //= BYTESIZE
            else: nBytes = nBytes//BYTESIZE + 1
            self.fp.write(pack("<L",int(nBytes))) # stores size as a little-endian unsigned long

            # create a PackedBits object to hold the nBytes of data for this channel/block of coded data
            pb = PackedBits()
            pb.Size(nBytes)

            # Write block type information (2 bits)
            if codingParams.currBlockType == LONG:
                pb.WriteBits(0, 2)  # 00 = LONG
            elif codingParams.currBlockType == START:
                pb.WriteBits(2, 2)  # 10 = START
            elif codingParams.currBlockType == SHORT:
                pb.WriteBits(1, 2)  # 01 = SHORT
            elif codingParams.currBlockType == STOP:
                pb.WriteBits(3, 2)  # 11 = STOP

            if codingParams.currBlockType == SHORT:
                num_short_blocks = 15

                for i in range(num_short_blocks):
                    # now pack the nBytes of data into the PackedBits object
                    pb.WriteBits(overallScaleFactor[iCh][i],codingParams.nScaleBits)  # overall scale factor
                    iMant=0  # index offset in mantissa array (because mantissas w/ zero bits are omitted)
                    for iBand in range(codingParams.sfBands.nBands): # loop over each scale factor band to pack its data
                        ba = bitAlloc[iCh][i][iBand]
                        if ba: ba-=1  # if non-zero, store as one less (since no bit allocation of 1 bits/mantissa)
                        pb.WriteBits(ba,codingParams.nMantSizeBits)  # bit allocation for this band (written as one less if non-zero)
                        pb.WriteBits(scaleFactor[iCh][i][iBand],codingParams.nScaleBits)  # scale factor for this band (if bit allocation non-zero)
                        if bitAlloc[iCh][i][iBand]:
                            for j in range(codingParams.sfBands.nLines[iBand]):
                                pb.WriteBits(mantissa[iCh][i][iMant+j],bitAlloc[iCh][i][iBand])     # mantissas for this band (if bit allocation non-zero) and bit alloc <>1 so is 1 higher than the number
                            iMant += codingParams.sfBands.nLines[iBand]  # add to mantissa offset if we passed mantissas for this band
            else:
                # now pack the nBytes of data into the PackedBits object
                pb.WriteBits(overallScaleFactor[iCh],codingParams.nScaleBits)  # overall scale factor
                iMant=0  # index offset in mantissa array (because mantissas w/ zero bits are omitted)
                for iBand in range(codingParams.sfBands.nBands): # loop over each scale factor band to pack its data
                    ba = bitAlloc[iCh][iBand]
                    if ba: ba-=1  # if non-zero, store as one less (since no bit allocation of 1 bits/mantissa)
                    pb.WriteBits(ba,codingParams.nMantSizeBits)  # bit allocation for this band (written as one less if non-zero)
                    pb.WriteBits(scaleFactor[iCh][iBand],codingParams.nScaleBits)  # scale factor for this band (if bit allocation non-zero)
                    if bitAlloc[iCh][iBand]:
                        for j in range(codingParams.sfBands.nLines[iBand]):
                            pb.WriteBits(mantissa[iCh][iMant+j],bitAlloc[iCh][iBand])     # mantissas for this band (if bit allocation non-zero) and bit alloc <>1 so is 1 higher than the number
                        iMant += codingParams.sfBands.nLines[iBand]  # add to mantissa offset if we passed mantissas for this band
            # done packing (end loop over scale factor bands)

            # CUSTOM DATA:
            # < now can add in custom data if space allocated in nBytes above>  # TODO for specific instructions on block? 
            # codingParams.secretMessage = "DATA BLOCK WRITING. BYE"
            # binary for block type
            # 0b00 = regular long, 0b01 = short, 0b10 = start, 0b11 = stop
            #codingParams.currblockType = bin(0)
            #codingParams.shortBlockInd = bin(7)  # indexed 0-7
            # finally, write the data in this channel's PackedBits object to the output file
            self.fp.write(pb.GetPackedData())
        # end loop over channels, done writing coded data for all channels
        return

    def Close(self,codingParams):
        """
        Flushes the last data block through the encoding process (if encoding)
        and closes the audio file
        """
        # determine if encoding or encoding and, if encoding, do last block
        if self.fp.mode == "wb":  # we are writing to the PACFile, must be encode
            # we are writing the coded file -- pass a block of zeros to move last data block to other side of MDCT block
            data = [ np.zeros(codingParams.nMDCTLines,dtype=np.float64),
                     np.zeros(codingParams.nMDCTLines,dtype=np.float64) ]
            self.WriteDataBlock(data, codingParams)
        self.fp.close()


    def Encode(self,data,codingParams):
        """
        Encodes multichannel audio data and returns a tuple containing
        the scale factors, mantissa bit allocations, quantized mantissas,
        and the overall scale factor for each channel.
        """
        #Passes encoding logic to the Encode function defined in the codec module
        return codec.Encode(data,codingParams)

    def Decode(self,scaleFactor,bitAlloc,mantissa, overallScaleFactor,codingParams):
        """
        Decodes a single audio channel of data based on the values of its scale factors,
        bit allocations, quantized mantissas, and overall scale factor.
        """
        #Passes decoding logic to the Decode function defined in the codec module
        return codec.Decode(scaleFactor,bitAlloc,mantissa, overallScaleFactor,codingParams)








#-----------------------------------------------------------------------------

# Testing the full PAC coder (needs a file called "input.wav" in the code directory)
if __name__=="__main__":

    print( "\nTesting the PAC coder (input.wav -> coded.pac -> output.wav):")
    import time
    from pcmfile import * # to get access to WAV file handling
    elapsed = time.time()

    for Direction in ("Encode", "Decode"):
#    for Direction in ("Decode",):

        # create the audio file objects
        if Direction == "Encode":
            print( "\n\tEncoding input PCM file...",)
            inFile= PCMFile("audio/spgm.wav")
            outFile = PACFile("audio/spgm_192kbps.pac")
        else: # "Decode"
            print( "\n\tDecoding coded PAC file...",)
            inFile = PACFile("audio/spgm_192kbps.pac")
            outFile= PCMFile("audio/spgm_192kbps.wav")
        # only difference is file names and type of AudioFile object

        # open input file
        codingParams=inFile.OpenForReading()  # (includes reading header)

        # pass parameters to the output file
        if Direction == "Encode":
            # set additional parameters that are needed for PAC file
            # (beyond those set by the PCM file on open)
            # codingParams.secretMessage = "starting encode"
            codingParams.nMDCTLines = 1024                            
            codingParams.nScaleBits = 3
            codingParams.nMantSizeBits = 5
            codingParams.targetBitsPerSample = 2.9
            # tell the PCM file how large the block size is
            codingParams.nSamplesPerBlock = codingParams.nMDCTLines   
        else: # "Decode"
            # set PCM parameters (the rest is same as set by PAC file on open)
            # codingParams.secretMessage = "starting decode"
            codingParams.bitsPerSample = 16
        # only difference is in setting up the output file parameters


        # open the output file
        outFile.OpenForWriting(codingParams) # (includes writing header)

        # Read the input file and pass its data to the output file to be written
        ''''
        1. prepend N zeros infront of first block, enter while loop. -- already included in starter code 
        2. readdatablock
        3. writedatablock:
            1. we will detect for transience in first block (current block). If find transience, 
            block type of zeros is start block. current block type is short. 
            2. do mdct on the zeros + first half of first block.
        - next iteration of while loop. repeat #1,2,3
            1. we will detect for transience in second block.
            2. do mdct on N/2 zeros + first half of first block.
        - next iteration of while loop. repeat #1,2,3, 3.1, 3.2
        '''
        # prepend N zeros infront of first block
        # data=[]
        # for iCh in range(codingParams.nChannels):
        #     data.append(np.zeros(codingParams.nMDCTLines, dtype=np.float64))
        # if true, read prepend zeros before start of actual audio data

        # used for detecting transience in audio file
        detector = TransientDetector(
            method=DetectionMethod.MODEL,
            model_path="transient_detection_model_final.pth"
        )
        # codingParams.firstBlock = True
        count = 0
        while True:
            # print("inside while loop: ", codingParams.__dict__)
            data = inFile.ReadDataBlock(codingParams)   # for encoding, see pcmfile ReadDataBlock
            # print("after read" , len(data[0]), " block")
            # print(codingParams.__dict__, "\n")
            if not data: break  # we hit the end of the input file
            outFile.WriteDataBlock(data,codingParams) 
            # if (codingParams.firstBlock): codingParams.firstBlock = False
            print( ".",end="")  # just to signal how far we've gotten to user
        # end loop over reading/writing the blocks

        # TODO: take out the prepended zeros to get back to original audio length
        # close the files
        inFile.Close(codingParams)
        outFile.Close(codingParams)
    # end of loop over Encode/Decode

    elapsed = time.time()-elapsed
    print( "\nDone with Encode/Decode test\n")
    print( elapsed ," seconds elapsed")
