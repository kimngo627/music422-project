"""

Music 422  Marina Bosi

window.py -- Defines functions to window an array of discrete-time data samples

-----------------------------------------------------------------------
© 2009-2025 Marina Bosi & Richard E. Goldberg -- All rights reserved
-----------------------------------------------------------------------


"""

### ADD YOUR CODE AT THE SPECIFIED LOCATIONS ###

import numpy as np
import scipy.fft
import matplotlib.pyplot as plt
import mdct

### Problem 1.d ###
def SineWindow(dataSampleArray):
    """
    Returns a copy of the dataSampleArray sine-windowed
    Sine window is defined following pp. 106-107 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    """

    ### YOUR CODE STARTS HERE ###
    N = len(dataSampleArray)
    window = np.sin(np.pi * (np.arange(N) + 0.5) / N)
    return dataSampleArray * window
    ### YOUR CODE ENDS HERE ###


def HanningWindow(dataSampleArray):
    """
    Returns a copy of the dataSampleArray Hanning-windowed
    Hann window is defined following pp. 106-107 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    """

    ### YOUR CODE STARTS HERE ###
    N = len(dataSampleArray)
    window = 0.5 * (1 - np.cos(2 * np.pi * (np.arange(N) + 0.5) / N))
    return dataSampleArray * window
    ### YOUR CODE ENDS HERE ###


### Problem 1.d - OPTIONAL ###
def KBDWindow(dataSampleArray,alpha=4.):
    """
    Returns a copy of the dataSampleArray KBD-windowed
    KBD window is defined following the KDB Window handout in the 
	Canvas Files/Assignments/HW3 folder
    """

    ### YOUR CODE STARTS HERE ###

    return np.zeros_like(dataSampleArray) # CHANGE THIS
    ### YOUR CODE ENDS HERE ###

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    ### YOUR TESTING CODE STARTS HERE ###
    N = 1024
    x = np.cos(2 * np.pi * 7000 * np.arange(N) / 44100)

    x_sine = SineWindow(x)
    fft_sine = np.abs(scipy.fft.fft(x_sine))
    mdct_sine = np.abs(mdct.MDCT(x_sine, N//2, N//2, isInverse=False))

    x_hanning = HanningWindow(x)
    fft_hanning = np.abs(scipy.fft.fft(x_hanning))

    fft_sine_dbspl = 96 + 10 * np.log10((4 / (N ** 2 * np.mean(x_sine ** 2))) * fft_sine ** 2)
    fft_hanning_dbspl = 96 + 10 * np.log10((4 / (N ** 2 * np.mean(x_hanning ** 2))) * fft_hanning ** 2)
    mdct_sine_dbspl = 96 + 10 * np.log10((2 / (np.mean(x_sine ** 2))) * (mdct_sine) ** 2)

    plt.figure(figsize=(12, 10))
    plt.subplot(2, 1, 1)
    plt.plot(x_sine, label="Sine Window")
    plt.plot(x_hanning, label="Hanning Window")
    plt.title("Windows in Time Domain")
    plt.legend()

    plt.subplot(2, 1, 2)
    plt.plot(fft_sine_dbspl, label="Sine Window (FFT)")
    plt.plot(fft_hanning_dbspl, label="Hanning Window (FFT)")
    plt.plot(mdct_sine_dbspl, label="Sine Window (MDCT)")
    plt.title("Frequency Content in dB SPL")
    plt.xlabel("Frequency Bin")
    plt.ylabel("dB SPL")
    plt.legend()

    plt.savefig('1f.png') 
    ### YOUR TESTING CODE ENDS HERE ###

