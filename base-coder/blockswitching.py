"""
blockswitching.py -- Defines functions to do window blockswitching on
an array of discrete-time data samples
* Using some concepts from Prof. Bosi & Goldberg's book. 
-----------------------------------------------------------------------
© 2009-2025 Marina Bosi & Richard E. Goldberg -- All rights reserved
-----------------------------------------------------------------------

"""
import numpy as np
import scipy.fft
import matplotlib.pyplot as plt
import mdct

def SineWindowFunc(N):
    '''N is length of the window'''
    return np.sin(np.pi * (np.arange(N) + 0.5) / N)

def HanningWindowFunc(N):
    '''N is length of the window'''
    return 0.5 * (1 - np.cos(2 * np.pi * (np.arange(N) + 0.5) / N))

def TransitionWindow(dataSampleArray, window, longLen, shortLen, start=True, short=False):
    '''
    See pg 123 in the book.
    Returns a copy of the dataSampleArray transition-windowed.
    Depending on if the window is the start or stop transition window,
    the shape of the window changes. 
    start window: [half of long, half of short]
    stop window: [half of shirt, half of long]
    @param dataSampleArray the data we apply windowing on
    @param window the type of windowing function we want to apply
    @param halfshortLen the length of half of a short block 
    @param halfLongLen the length of half of a long block
    @param start determines what kind of transition window is applied
    '''
    N = len(dataSampleArray) 
    halfLongLen = longLen//2
    halfShortLen = shortLen//2
    # print(f"start? {start} N: {N}, halfLongLen: {halfLongLen}, halfShortLen: {halfShortLen}")
    if (N < halfLongLen):
        print("not enough data samples to window with half LONG window")
    elif (N < halfShortLen):
        print("not enough data samples to window with half SHORT window")
    if (short):
        return window(shortLen) * dataSampleArray[:shortLen]
    if (start):   #half long on left, half short on right
        halfLongWindowed = window(longLen)[:halfLongLen] * dataSampleArray[:halfLongLen]
        # Create a smooth transition window between long and short
        transition_len = halfLongLen
        transition_window = np.zeros(transition_len)
        # This uses a cos^2 + sin^2 = 1 identity for perfect reconstruction
        for n in range(transition_len):
            alpha = n / transition_len
            transition_window[n] = np.sqrt(1 - (window(longLen)[halfLongLen + n] ** 2)) 
        
        halfShortWindowed = window(shortLen)[halfShortLen:] * dataSampleArray[halfLongLen:halfLongLen+halfShortLen]
        
        return np.concatenate((halfLongWindowed, halfShortWindowed))
    else:  #half short on left, half long on right
        halfShortWindowed = window(shortLen)[:halfShortLen] * dataSampleArray[:halfShortLen]
        halfLongWindowed = window(longLen)[halfLongLen:] * dataSampleArray[halfShortLen:halfShortLen+halfLongLen]
        
        return np.concatenate((halfShortWindowed, halfLongWindowed))



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

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":
    N = 1024
    x = np.cos(2 * np.pi * 7000 * np.arange(N) / 44100)
    shortN = N//2
    # short_x_firsthalf = np.cos(2 * np.pi * 7000 * np.arange(shortN) / 44100)
    # short_x_sechalf = np.cos(2 * np.pi * 7000 * np.arange(shortN, shortN+N//2) / 44100)
    short_x_half = [1]*(N//2)
    ones = [1]*N
    ############################ Testing Transition Window Implementation ############################
    # usually don't use hanning window for transition window, use it for fft
    # x_hann = HanningWindow(ones)
    # short_x_hann = np.concat((HanningWindow(short_x_firsthalf), HanningWindow(short_x_sechalf)))
    # x_hann_transition_start = TransitionWindow(ones, HanningWindowFunc, N//2, N//4, True)
    # x_hann_transition_stop = TransitionWindow(ones, HanningWindowFunc, N//2, N//4, False)
    # plt.plot(x_hann, label="Long Hann Window")
    # plt.plot(short_x_hann, label="Short Hann Window")
    # plt.plot(x_hann_transition_start, label="Start Transition Hanning Window Ones")
    # plt.plot(x_hann_transition_stop, label="Stop Transition Hanning Window")

    x_sine = SineWindow(ones)
    short_x_sine = np.concat((SineWindow(short_x_half) , SineWindow(short_x_half)))
    x_sine_transition_start = TransitionWindow(ones, SineWindowFunc, N//2, N//4, True)
    x_sine_transition_stop = TransitionWindow(ones, SineWindowFunc, N//2, N//4, False)
    
    plt.figure(figsize=(12, 10))
    plt.subplot(2, 1, 1)
    plt.title("Regular Long Windows")
    plt.plot(x_sine, label="Long Sine Window")
    
    plt.subplot(2, 1, 2)
    plt.title("Regular Short Windows")
    plt.plot(short_x_sine, label="Short Sine Window")
    
    plt.xlabel('Time (samples)')
    plt.ylabel('Magnitude')
    plt.legend()
    plt.savefig('Long vs Short Window')
    plt.show()

    plt.figure(figsize=(12, 10))
    plt.subplot(3, 1, 1)
    plt.title("Regular Long Windows")
    plt.plot(x_sine, label="Regular Sine Window")
    
    plt.subplot(3, 1, 2)
    plt.plot(x_sine_transition_start, label="Start Transition Sine Window")
    plt.title("Start Window in Time Domain")
    plt.axvline(x=N//4, color='grey', linestyle='--', label=f'{N//4}')
    plt.legend()

    plt.subplot(3, 1, 3)
    plt.plot(x_sine_transition_stop, label="Stop Transition Sine Window")
    plt.title('Stop Window in Time Domain')
    plt.axvline(x=N//4//2, color='grey', linestyle='--', label=f'{N//8}')
    plt.xlabel('Time (samples)')
    plt.ylabel('Magnitude')
    plt.legend()
    plt.subplots_adjust(wspace=0.4, hspace=0.4)  # Adjust values as needed
    plt.savefig('Transition Windows Ones') 
    plt.show()
    ####################################### windowing on signal x #######################################
    # x_sine = SineWindow(x)
    # x_hann = HanningWindow(x)
    # short_x_sine = np.concat((SineWindow(short_x_firsthalf) , SineWindow(short_x_sechalf)))
    # short_x_hann = np.concat((HanningWindow(short_x_firsthalf), HanningWindow(short_x_sechalf)))
    # x_sine_transition_start = TransitionWindow(x, SineWindowFunc, N//2, N//4, True)
    # x_hann_transition_start = TransitionWindow(x, HanningWindowFunc, N//2, N//4, True)

    # x_sine_transition_stop = TransitionWindow(x, SineWindowFunc, N//2, N//4, False)
    # x_hann_transition_stop = TransitionWindow(x, HanningWindowFunc, N//2, N//4, False)

    # plt.figure(figsize=(12, 10))
    # plt.subplot(2, 1, 1)
    # plt.title("Regular Long Windows")
    # plt.plot(x_sine, label="Long Sine Window")
    # plt.plot(x_hann, label="Long Hann Window")
    # plt.subplot(2, 1, 2)
    # plt.title("Regular Short Windows")
    # plt.plot(short_x_sine, label="Short Sine Window")
    # plt.plot(short_x_hann, label="Short Hann Window")
    # plt.xlabel('Time (samples)')
    # plt.ylabel('Magnitude')
    # plt.legend()
    # plt.savefig('Long vs Short Window')
    # plt.show()

    # plt.figure(figsize=(12, 10))
    # plt.subplot(3, 1, 1)
    # plt.title("Regular Long Windows")
    # plt.plot(x_sine, label="Regular Sine Window")
    # plt.plot(x_hann, label="Regular Hann Window")
    # plt.subplot(3, 1, 2)
    # plt.plot(x_sine_transition_start, label="Start Transition Sine Window")
    # plt.plot(x_hann_transition_start, label="Start Transition Hanning Window")
    # plt.title("Start Windows in Time Domain")
    # plt.legend()

    # plt.subplot(3, 1, 3)
    # plt.plot(x_sine_transition_stop, label="Stop Transition Sine Window")
    # plt.plot(x_hann_transition_stop, label="Stop Transition Hanning Window")
    # plt.title('Stop Windows in Time Domain')
    # plt.xlabel('Time (samples)')
    # plt.ylabel('Magnitude')
    # plt.legend()
    # plt.subplots_adjust(wspace=0.4, hspace=0.4)  # Adjust values as needed
    # plt.savefig('Transition Windows') 
    # plt.show()
    ############################ PLOTTING WINDOWS ####################################################
    # x_sine = SineWindow(x)
    # fft_sine = np.abs(scipy.fft.fft(x_sine))
    # mdct_sine = np.abs(mdct.MDCT(x_sine, N//2, N//2, isInverse=False))

    # x_hanning = HanningWindow(x)
    # fft_hanning = np.abs(scipy.fft.fft(x_hanning))

    # fft_sine_dbspl = 96 + 10 * np.log10((4 / (N ** 2 * np.mean(x_sine ** 2))) * fft_sine ** 2)
    # fft_hanning_dbspl = 96 + 10 * np.log10((4 / (N ** 2 * np.mean(x_hanning ** 2))) * fft_hanning ** 2)
    # mdct_sine_dbspl = 96 + 10 * np.log10((2 / (np.mean(x_sine ** 2))) * (mdct_sine) ** 2)

    # plt.figure(figsize=(12, 10))
    # plt.subplot(2, 1, 1)
    # plt.plot(x_sine, label="Sine Window")
    # plt.plot(x_hanning, label="Hanning Window")
    # plt.title("Windows in Time Domain")
    # plt.legend()

    # plt.subplot(2, 1, 2)
    # plt.plot(fft_sine_dbspl, label="Sine Window (FFT)")
    # plt.plot(fft_hanning_dbspl, label="Hanning Window (FFT)")
    # plt.plot(mdct_sine_dbspl, label="Sine Window (MDCT)")
    # plt.title("Frequency Content in dB SPL")
    # plt.xlabel("Frequency Bin")
    # plt.ylabel("dB SPL")
    # plt.legend()

    # plt.savefig('1f.png') 
    ################################################################################################

