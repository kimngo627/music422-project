"""
Music 422
-----------------------------------------------------------------------
(c) 2009-2025 Marina Bosi  -- All rights reserved
-----------------------------------------------------------------------
"""

import numpy as np
from window import *
import matplotlib.pyplot as plt
from psychoac import *
from mdct import *
from quantize import * 
import scipy.fft
import scipy.signal

# Question 1.c)
def BitAllocUniform(bitBudget, maxMantBits, nBands, nLines, SMR=None):
    """
    Returns a hard-coded vector that, in the case of the signal used in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are uniformly distributed for the mantissas.
    """
    N = np.sum(nLines)
    bits = np.zeros(nBands)
    
    while(bitBudget >= np.min(nLines)):
        for k in range(nBands):
            if(nLines[k] < bitBudget):
                bits[k] += 1
                bitBudget -= nLines[k]
        
    return bits.astype(np.int64)
    #return np.array([3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2])

def BitAllocConstSNR(bitBudget, maxMantBits, nBands, nLines, peakSPL):
    """
    Returns a hard-coded vector that, in the case of the signal used in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are distributed for the mantissas to try and keep a constant
    quantization noise floor (assuming a noise floor 6 dB per bit below
    the peak SPL line in the scale factor band).
    """
    bits = np.zeros(nBands)
    SNR = peakSPL
    i = 0
    
    while(bitBudget >= np.min(nLines)):
        i = np.argmax(SNR)
        if(bitBudget >= nLines[i]):
            bits[i] += 1
            bitBudget -= nLines[i]
        SNR[i] -= 6.0
        
    single_bits = np.where(bits == 1)
    bits[single_bits] -= 1
    bitBudget += np.sum(nLines[single_bits])
    
    found_max = np.where(bits > maxMantBits)
    extra_bits = bits[found_max] - maxMantBits
    bits[found_max] -= extra_bits
    bitBudget += np.sum(nLines[found_max]*extra_bits)
        
    while(bitBudget >= np.min(nLines)):
        i = np.argmax(SNR)
        if(bitBudget >= nLines[i] and bits[i] > 0 and bits[i] < maxMantBits):
            bits[i] += 1
            bitBudget -= nLines[i]
        SNR[i] -= 6.0
            
    return bits.astype(np.int64) 
    return np.array([9, 11, 16, 13, 14, 10, 8, 13, 9, 6, 5, 4, 4, 3, 3, 2, 2, 8, 11, 0, 0, 10, 0, 0, 0])

def BitAllocConstNMR(bitBudget, maxMantBits, nBands, nLines, SMR):
    """
    Returns a hard-coded vector that, in the case of the signal used in HW#4,
    gives the allocation of mantissa bits in each scale factor band when
    bits are distributed for the mantissas to try and keep the quantization
    noise floor a constant distance below (or above, if bit starved) the
    masked threshold curve (assuming a quantization noise floor 6 dB per
    bit below the peak SPL line in the scale factor band).
    """
    bits = np.zeros(nBands)
    
    while(bitBudget >= np.min(nLines)):    
        i = np.argmax(SMR)
        if(bitBudget >= nLines[i]):
            bits[i] += 1
            bitBudget -= nLines[i]
        SMR[i] -= 6.0
        
    single_bit_indices = np.where(bits == 1)
    bits[single_bit_indices] -= 1
    bitBudget += np.sum(nLines[single_bit_indices])

    found_max = np.where(bits > maxMantBits)
    extra_bits = bits[found_max] - maxMantBits
    bits[found_max] -= extra_bits
    bitBudget += np.sum(nLines[found_max]*extra_bits)
        
    while(bitBudget >= np.min(nLines)):
        i = np.argmax(SMR)
        if(bitBudget >= nLines[i] and bits[i] > 0 and bits[i] < maxMantBits):
            bits[i] += 1
            bitBudget -= nLines[i]
        SMR[i] -= 6.0

    return bits.astype(np.int64)
    return np.array([7, 7, 9, 9, 10, 7, 7, 9, 5, 6, 7, 9, 9, 9, 8, 8, 4, 4, 8, 0, 2, 8, 0, 0, 0])

# Question 2.a)
def BitAlloc(bitBudget, maxMantBits, nBands, nLines, SMR):
    """
    Allocates bits to scale factor bands so as to flatten the NMR across the spectrum

       Arguments:
           bitBudget is total number of mantissa bits to allocate
           maxMantBits is max mantissa bits that can be allocated per line
           nBands is total number of scale factor bands
           nLines[nBands] is number of lines in each scale factor band
           SMR[nBands] is signal-to-mask ratio in each scale factor band

        Returns:
            bits[nBands] is number of bits allocated to each scale factor band

        
    """
    priorities = np.log2(np.clip(SMR, 1e-10, None))
    sorted_indices = np.argsort(priorities)[::-1]
    bits = np.zeros(nBands, dtype=int)

    threshold = np.max(priorities)

    while bitBudget >= np.min(nLines):
        allocated = False 

        for i in sorted_indices:
            if priorities[i] >= threshold and bits[i] < maxMantBits and bitBudget >= nLines[i]:
                bits[i] += 1
                bitBudget -= nLines[i]
                allocated = True

                if bitBudget <= 0:
                    break

        if not allocated:
            break 

        threshold -= 1

    single_bit_indices = np.where(bits == 1)[0]
    bitBudget += np.sum(nLines[single_bit_indices])
    bits[single_bit_indices] = 0 

    for i in sorted_indices:
        if bits[i] < maxMantBits and bitBudget >= nLines[i]:
            if bits[i] == 0 and bitBudget >= 2 * nLines[i]:
                bits[i] += 2
                bitBudget -= 2 * nLines[i]
                allocated = True
            elif bits[i] == 0 and bitBudget < 2 * nLines[i]:
                allocated = True
            else:
                bits[i] += 1
                bitBudget -= nLines[i]
                allocated = True

            if bitBudget <= 0:
                break 

    return bits.astype(np.int64)


#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    Fs = 48000
    amp = np.array([0.40, 0.24, 0.18, 0.08, 0.04, 0.02])
    freqs = np.array([220, 330, 440, 880, 4400, 8800], dtype = float)
    L = 1024
    x = np.zeros(L, dtype = float)
    n = np.arange(L)    
    
    for k in range(6):
        x = x + amp[k]*np.cos(2*np.pi*freqs[k]*n/Fs)
         
    N = 1024
    x_win = HanningWindow(x)
    X_fft = scipy.fft.fft(x_win)
    X_fft = X_fft[:N//2]
    data = SineWindow(x)
    MDCTdata = MDCT(data, N//2, N//2)
    MDCT_spl = np.maximum(96. + 10*np.log10(4.*MDCTdata**2),-30)
    freqs_mdct = np.linspace(0,Fs/2,N//2+1) + 0.5*(Fs/N)

    freqs_mdct = freqs_mdct[:N//2]
               
    peak_SPL, _ = scipy.signal.find_peaks(MDCT_spl)
    
    
    nBands = 25
    dataRate = 192000
    nMDCTLines = N/2.0
    bitBudget = (dataRate * nMDCTLines)/Fs - 2*(nBands*4) - 4
    print(bitBudget)
    sfBands = ScaleFactorBands(AssignMDCTLinesFromFreqLimits(nMDCTLines,Fs))

    MDCTscale = 1 
    masked_threshold = getMaskedThreshold(x, MDCTdata, MDCTscale, Fs, sfBands)
    mdct_spl = SPL(4 * MDCTdata**2)
    
    bits_uniform = BitAllocUniform(bitBudget, 16, nBands, sfBands.nLines)
    
    peakSPL = np.zeros(nBands)
    for k in range(nBands):
        peakSPL[k] = np.max(MDCT_spl[sfBands.lowerLine[k]:sfBands.upperLine[k]+1])
 
    bits_constSNR = BitAllocConstSNR(bitBudget, 16, nBands, sfBands.nLines, peakSPL.copy())
    
    SMR = CalcSMRs(data, MDCTdata, 0, Fs, sfBands)

    bits_constNMR = BitAllocConstNMR(bitBudget, 16, nBands, sfBands.nLines, SMR.copy())

    bits_opt = BitAlloc(bitBudget, 16, nBands, sfBands.nLines, SMR.copy())
    
    zVec = Bark(freqs_mdct)    
    mask = list()
    intensity = list(np.array([]))
    
    for i in range(len(peak_SPL)):
        mask.append(Masker(freqs_mdct[i], peak_SPL[i], True))
        intensity.append(mask[i].vIntensityAtBark(zVec))
        
        
    bits = [bits_uniform, bits_constSNR, bits_constNMR, bits_opt]
    freqs_bark = np.array([100,200,300,400,510,630,770,920,1080,1270,1400,1720,2000,2320,2700,3150, 3700,4400,5300,6400,7700,9500,12000,15500, 24000])
    cb_barks = Bark(freqs_bark)
    strategy = ["Uniform", "SNR", "NMR", "Water-Filling"]
    for i in range(4):
        noise_floor = peakSPL - 6.0 * bits[i]
        plt.figure(i+1, figsize=(10, 5))
        for cb_z in freqs_bark:
            plt.axvline(x=cb_z, color="gray", linestyle="dashed", alpha=0.5)
        plt.step(freqs_mdct, MDCT_spl, label="MDCT SPL")
        plt.step(freqs_mdct, masked_threshold, 'k', label="Masked Threshold")
        plt.step(freqs_bark, noise_floor, label=f"{strategy[i]} Noise Floor")
        plt.ylim(ymax = 100)
        plt.xlim(xmin = 10, xmax = 16000)
        plt.xlabel('Log frequency in Hz')
        plt.ylabel('Normalized SPL in dB')
        plt.legend()
        plt.grid()
        plt.show()

    print("\nCritical Band Bit Allocation:")
    print(f"{'Band':<6}{'Uniform':<12}{'ConstSNR':<15}{'ConstNMR':<15}{'Water-Filling':<15}")
    for i in range(sfBands.nBands):
        print(f"{i:<6}{bits_uniform[i]:<12}{bits_constSNR[i]:<15}{bits_constNMR[i]:<15}{bits_opt[i]:<15}")
