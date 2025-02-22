"""
psychoac.py -- masking models implementation

-----------------------------------------------------------------------
(c) 2011-2025 Marina Bosi & Richard E. Goldberg -- All rights reserved
-----------------------------------------------------------------------
"""

import numpy as np
from window import *
import scipy.fft
import scipy.signal
from mdct import *

def SPL(intensity):
    """
    Returns the SPL corresponding to intensity 
    """
    intensity = np.maximum(intensity, 1e-100)
    spl = 96 + 10 * np.log10(intensity)
    return np.maximum(spl, -30)

def Intensity(spl):
    """
    Returns the intensity  for SPL spl
    """
    return 10 ** ((spl - 96) / 10)

def Thresh(f):
    """Returns the threshold in quiet measured in SPL at frequency f (in Hz)"""

    f = np.maximum(f, 20)

    thresh = 3.64 * (f / 1000)**(-0.8) - 6.5 * np.exp(-0.6 * ((f / 1000) - 3.3)**2) + 10**(-3) * (f / 1000) ** 4
    return thresh

def Bark(f):
    """Returns the bark-scale frequency for input frequency f (in Hz) """
    return 13 * np.arctan(0.76 * f / 1000) + 3.5 * np.arctan((f / 7500) ** 2)

class Masker:
    """
    a Masker whose masking curve drops linearly in Bark beyond 0.5 Bark from the
    masker frequency
    """

    def __init__(self,f,SPL,isTonal=True):
        """
        initialized with the frequency and SPL of a masker and whether or not
        it is Tonal
        """
        self.f = f
        self.SPL = SPL
        self.delta = 16 if isTonal else 6
        self.z = Bark(f)

    def IntensityAtFreq(self,freq):
        """The intensity at frequency freq"""
        return 0 # TO REPLACE WITH YOUR CODE

    def IntensityAtBark(self,z):
        """The intensity at Bark location z"""
        dz = z - self.z

        if abs(dz) <= 0.5:
            spread = 0
        elif dz < -0.5:
            spread = -27 * (abs(dz) - 0.5)
        else:
            spread = (-27 + 0.367 * max(self.SPL - 40, 0)) * (abs(dz) - 0.5)
        
        return self.SPL - self.delta + spread

    def vIntensityAtBark(self,zVec):
        """The intensity at vector of Bark locations zVec"""

        dz = zVec - self.z
        spread = np.where(
            np.abs(dz) <= 0.5, 0,
            np.where(dz < -0.5, -27 * (np.abs(dz) - 0.5),
                     (-27 + 0.367 * max(self.SPL - 40, 0)) * (np.abs(dz) - 0.5))
        )
        return self.SPL - self.delta + spread

        return intensities


# Default data for 25 scale factor bands based on the traditional 25 critical bands
cbFreqLimits = np.array([100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000,
    2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000, 15500, 24000])  # TO REPLACE WITH THE APPROPRIATE VALUES

def AssignMDCTLinesFromFreqLimits(nMDCTLines, sampleRate, flimit = cbFreqLimits):
    """
    Assigns MDCT lines to scale factor bands for given sample rate and number
    of MDCT lines using predefined frequency band cutoffs passed as an array
    in flimit (in units of Hz). If flimit isn't passed it uses the traditional
    25 Zwicker & Fastl critical bands as scale factor bands.
    """
    bin_width = (sampleRate / 2) / nMDCTLines
    maxFreq = (nMDCTLines - 1 + 0.5) * bin_width
    nLines = []
    iLast = -1
    for iLine in range(len(flimit)):
        if flimit[iLine] > maxFreq:
            nLines.append(nMDCTLines - 1 - iLast)
            break
        iUpper = int(flimit[iLine] / bin_width - 0.5)
        nLines.append(iUpper - iLast)
        iLast = iUpper
    return nLines

class ScaleFactorBands:
    """
    A set of scale factor bands (each of which will share a scale factor and a
    mantissa bit allocation) and associated MDCT line mappings.

    Instances know the number of bands nBands; the upper and lower limits for
    each band lowerLimit[i in range(nBands)], upperLimit[i in range(nBands)];
    and the number of lines in each band nLines[i in range(nBands)]
    """

    def __init__(self,nLines):
        """
        Assigns MDCT lines to scale factor bands based on a vector of the number
        of lines in each band
        """
        self.nBands = len(nLines)
        self.nLines = np.array(nLines, dtype=np.uint16)
        self.lowerLine = np.empty(self.nBands, dtype=np.uint16)
        self.upperLine = np.empty(self.nBands, dtype=np.uint16)
        self.lowerLine[0] = 0
        self.upperLine[0] = nLines[0] - 1
        for iBand in range(1, self.nBands):
            self.lowerLine[iBand] = self.upperLine[iBand - 1] + 1
            self.upperLine[iBand] = self.upperLine[iBand - 1] + nLines[iBand]



def getMaskedThreshold(data, MDCTdata, MDCTscale, sampleRate, sfBands):
    """
    Return Masked Threshold evaluated at MDCT lines.

    Used by CalcSMR, but can also be called from outside this module, which may
    be helpful when debugging the bit allocation code.
    """

    N = len(data)
    nLines = N//2
    lineToFreq = sampleRate / N
    fftData = scipy.fft.fft(HanningWindow(data))[:nLines]
    fftIntensity = 32./3./N/N*np.abs(fftData)**2

    maskers = []

    for i in range(2, nLines - 2):
        if fftIntensity[i]>fftIntensity[i-1] and fftIntensity[i]>fftIntensity[i+1]:
            spl = fftIntensity[i] + fftIntensity[i-1]+fftIntensity[i+1]

            f = lineToFreq * (i * fftIntensity[i] + (i-1) * fftIntensity[i-1] + (i+1) * fftIntensity[i+1]) / spl

            spl = SPL(spl)

            if spl>Thresh(f):
                maskers.append(Masker(f, spl))

    fline = lineToFreq*np.linspace(0.5, nLines - 0.5, nLines)
    zline = Bark(fline)

    maskedSPL = np.zeros(nLines, dtype=np.float64)

    for m in maskers: maskedSPL += Intensity(m.vIntensityAtBark(zline))
    maskedSPL += Intensity(Thresh(fline))
    maskedSPL = SPL(maskedSPL)
    return maskedSPL

def CalcSMRs(data, MDCTdata, MDCTscale, sampleRate, sfBands):
    """
    Set SMR for each critical band in sfBands.

    Arguments:
                data:       is an array of N time domain samples
                MDCTdata:   is an array of N/2 MDCT frequency coefficients for the time domain samples
                            in data; note that the MDCT coefficients have been scaled up by a factor
                            of 2^MDCTscale
                MDCTscale:  corresponds to an overall scale factor 2^MDCTscale for the set of MDCT
                            frequency coefficients
                sampleRate: is the sampling rate of the time domain samples
                sfBands:    points to information about which MDCT frequency lines
                            are in which scale factor band

    Returns:
                SMR[sfBands.nBands] is the maximum signal-to-mask ratio in each
                                    scale factor band

    Logic:
                Performs an FFT of data[N] and identifies tonal and noise maskers.
                Combines their relative masking curves and the hearing threshold
                to calculate the overall masked threshold at the MDCT frequency locations. 
				Then determines the maximum signal-to-mask ratio within
                each critical band and returns that result in the SMR[] array.
    """
    masked_threshold = getMaskedThreshold(data, MDCTdata, MDCTscale, sampleRate, sfBands)

    N_mdct = len(MDCTdata)
    freq_bins = np.linspace(0, sampleRate / 2, N_mdct)
    window_energy = 1/2

    MDCTdata_scaled = MDCTdata * (2 ** (-MDCTscale)) / (2/N_mdct)

    intensity_spectrum = (8 / (N_mdct**2 * window_energy)) * (np.abs(MDCTdata_scaled) ** 2)
    mdct_spl = SPL(intensity_spectrum)

    # plt.figure(figsize=(12, 6))
    # plt.plot(freq_bins, mdct_spl, label="MDCT SPL", linewidth=1)
    # plt.plot(freq_bins, masked_threshold, label="Masking Threshold", linestyle="dashed", color="red")
    # plt.xscale("log")
    # plt.xlabel("Frequency (Hz)")
    # plt.xlim(0, 10000)
    # plt.ylabel("SPL (dB)")
    # plt.ylim(-20, 100)
    # plt.title("MDCT Spectrum vs. Masking Threshold")
    # plt.legend()
    # plt.grid()
    # plt.show()

    smr_values = np.zeros(sfBands.nBands)

    for iBand in range(sfBands.nBands):
        band_start = sfBands.lowerLine[iBand]
        band_end = sfBands.upperLine[iBand] + 1

        smr_values[iBand] = np.max(mdct_spl[band_start:band_end] - masked_threshold[band_start:band_end])

    return smr_values

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":
    
    # 1.b
    # magnitude -> intensity -> spl
    Fs = 48000  
    frequencies = [220, 330, 440, 880, 4400, 8800]  
    amplitudes = [0.40, 0.24, 0.18, 0.08, 0.04, 0.02]  
    block_sizes = [512, 1024, 2048]  
    time_duration = 1 

    for N in block_sizes:
        t = np.arange(N) / Fs

        x = sum(A * np.cos(2 * np.pi * f * t) for A, f in zip(amplitudes, frequencies))

        windowed_signal = HanningWindow(x)
        window_energy = 3/8

        X = scipy.fft.fft(windowed_signal)
        magnitude_spectrum = np.abs(X[:N // 2]) ** 2 

        intensity_spectrum = (4 / (N**2 * window_energy)) * magnitude_spectrum

        spl_spectrum = SPL(intensity_spectrum)

        freq_bins = np.linspace(0, Fs / 2, N // 2)

        peak_indices, _ = scipy.signal.find_peaks(spl_spectrum)

        refined_frequencies = []
        refined_spl = []
        for kp in peak_indices:
            num = sum(m * intensity_spectrum[m] for m in range(kp - 1, kp + 2))
            denom = sum(intensity_spectrum[m] for m in range(kp - 1, kp + 2))
            fp = (Fs / N) * (num / denom)
            refined_frequencies.append(fp)
            num_spl = sum(spl_spectrum[m] for m in range(kp - 1, kp + 2))
            denom_spl = 3
            splp = (num_spl / denom_spl)
            refined_spl.append(splp)

        print(f"\nResults for N = {N}:")
        print("Peak Frequency (Hz) | Estimated SPL (dB)")
        #
        for f, spl_val in zip(refined_frequencies, refined_spl):
            print(f"{f:10.2f} Hz | {spl_val:10.2f} dB")

        plt.figure(figsize=(10, 5))
        plt.plot(freq_bins, spl_spectrum, label=f"N={N}", linewidth=1)
        plt.scatter(refined_frequencies, refined_spl, color='red', label="Estimated Peaks")
        plt.scatter(freq_bins[peak_indices], spl_spectrum[peak_indices], color='orange', label="Actual Peaks")
        plt.xscale("log")  
        plt.xlim(0, 10000)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("SPL (dB)")
        plt.title(f"SPL Spectrum for N={N}")
        plt.legend()
        plt.grid()
        plt.show()

    # 1.c
    Fs = 48000  
    N = 1024  
    frequencies = [220, 330, 440, 880, 4400, 8800]  
    amplitudes = [0.40, 0.24, 0.18, 0.08, 0.04, 0.02]  

    t = np.arange(N) / Fs
    test_signal = sum(A * np.cos(2 * np.pi * f * t) for A, f in zip(amplitudes, frequencies))

    windowed_signal = HanningWindow(test_signal)
    window_energy = 3/8

    X = scipy.fft.fft(windowed_signal)
    magnitude_spectrum = np.abs(X[:N // 2]) ** 2 
    intensity_spectrum = (4 / (N**2 * window_energy)) * magnitude_spectrum
    spl_spectrum = SPL(intensity_spectrum)

    freq_bins = np.linspace(0, Fs / 2, N // 2)

    quiet_threshold = Thresh(freq_bins)

    plt.figure(figsize=(10, 5))
    plt.plot(freq_bins, spl_spectrum, label="SPL of Test Signal", linewidth=1)
    plt.plot(freq_bins, quiet_threshold, label="Threshold in Quiet", linestyle="dashed", color="red")
    plt.xscale("log")  
    plt.xlabel("Frequency (Hz)")
    plt.xlim(0, 10000)
    plt.ylabel("SPL (dB)")
    plt.ylim(-20, 100)
    plt.title("SPL Spectrum vs. Threshold in Quiet")
    plt.legend()
    plt.grid()
    plt.show()

    # 1.d
    critical_band_frequencies = np.array([
        0, 100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720,
        2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000, 15500
    ])

    bark_values = Bark(critical_band_frequencies)

    print("Lower Frequency (Hz) | Computed Bark Value")
    for f, z in zip(critical_band_frequencies, bark_values):
        print(f"{f:18.1f} | {z:5.2f}")

    # 1.e
    Fs = 48000 
    N = 1024  
    frequencies = [220, 330, 440, 880, 4400, 8800]  
    amplitudes = [0.40, 0.24, 0.18, 0.08, 0.04, 0.02] 

    t = np.arange(N) / Fs
    test_signal = sum(A * np.cos(2 * np.pi * f * t) for A, f in zip(amplitudes, frequencies))


    windowed_signal = HanningWindow(test_signal)
    window_energy = 3/8


    X = scipy.fft.fft(windowed_signal)
    magnitude_spectrum = np.abs(X[:N // 2]) ** 2  
    intensity_spectrum = (4 / (N**2 * window_energy)) * magnitude_spectrum
    spl_spectrum = SPL(intensity_spectrum)

    freq_bins = np.linspace(0, Fs / 2, N // 2)

    quiet_threshold = Thresh(freq_bins)

    peak_indices, _ = scipy.signal.find_peaks(spl_spectrum)
    peak_frequencies = freq_bins[peak_indices]
    peak_spl_values = spl_spectrum[peak_indices]

    maskers = [Masker(f, spl) for f, spl in zip(peak_frequencies, peak_spl_values)]

    bark_scale = Bark(freq_bins)
    masking_curves = np.array([masker.vIntensityAtBark(bark_scale) for masker in maskers])
    masked_threshold = np.maximum.reduce(masking_curves) 
    masked_threshold = np.maximum(masked_threshold, quiet_threshold)

    plt.figure(figsize=(10, 6))

    bark_bins = Bark(freq_bins) 

    plt.plot(bark_bins, spl_spectrum, label="SPL of Test Signal", linewidth=1)

    plt.plot(bark_bins, quiet_threshold, label="Threshold in Quiet", linestyle="dotted", color="red")

    for masker, masking_curve in zip(maskers, masking_curves):
        plt.plot(bark_bins, masking_curve, linestyle="solid", label=f"Masker at {masker.f:.1f} Hz")

    plt.plot(bark_bins, masked_threshold, linestyle="solid", color="black", label="Masked Threshold")

    cb_freqs = np.array([
        100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000,
        2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000, 15500
    ])

    cb_barks = Bark(cb_freqs)

    for cb_z in cb_barks:
        plt.axvline(x=cb_z, color="gray", linestyle="dashed", alpha=0.5)

    plt.xlabel("Bark Scale")
    plt.xlim(0, 24) 
    plt.ylabel("SPL (dB)")
    plt.ylim(-20, 100)
    plt.title("Masking Thresholds and SPL in Bark Scale")
    plt.legend()
    plt.grid()
    plt.show()

    # 1.f
    nMDCTLines = 512 
    Fs = 48000 

    nLines = AssignMDCTLinesFromFreqLimits(nMDCTLines, Fs)

    sfBands = ScaleFactorBands(nLines)

    print("\nScale Factor Bands (N=1024, Fs=48kHz):")
    print(f"{'Band':<5} {'Lower Limit':<12} {'Upper Limit':<12} {'Lines':<8}")
    print("-" * 40)
    for i in range(sfBands.nBands):
        print(f"{i:<5} {sfBands.lowerLine[i]:<12} {sfBands.upperLine[i]:<12} {sfBands.nLines[i]:<8}")

    # 1.g
    Fs = 48000  
    N = 1024  
    MDCTscale = 1

    frequencies = [220, 330, 440, 880, 4400, 8800]  
    amplitudes = [0.40, 0.24, 0.18, 0.08, 0.04, 0.02] 

    t = np.arange(N) / Fs
    test_signal = sum(A * np.cos(2 * np.pi * f * t) for A, f in zip(amplitudes, frequencies))
    test_signal_sine = SineWindow(test_signal)

    MDCTdata = MDCT(test_signal_sine, N//2, N//2)[:N//2]
    MDCTdata_scaled = MDCTdata * (2 ** MDCTscale)

    sfBands = ScaleFactorBands(AssignMDCTLinesFromFreqLimits(N // 2, Fs))

    smr_values = CalcSMRs(test_signal, MDCTdata_scaled, MDCTscale, Fs, sfBands)

    bin_width = Fs / (2 * (N // 2))

    plt.figure(figsize=(12, 6))
    plt.bar(range(sfBands.nBands), smr_values, color='blue', alpha=0.7)
    plt.xlabel("Scale Factor Band Index")
    plt.ylabel("SMR (dB)")
    plt.title("Signal-to-Mask Ratio (SMR) for Each Scale Factor Band")
    plt.grid()
    plt.show()

    print("\nScale Factor Band SMR Table:")
    print(f"{'Band Index':<12}{'Lower Freq (Hz)':<18}{'Upper Freq (Hz)':<18}{'SMR (dB)':<10}")
    print("-" * 60)

    for iBand in range(sfBands.nBands):
        lower_freq = sfBands.lowerLine[iBand] * bin_width
        upper_freq = sfBands.upperLine[iBand] * bin_width
        smr_value = smr_values[iBand]
        print(f"{iBand:<12}{lower_freq:<18.2f}{upper_freq:<18.2f}{smr_value:<10.2f}")