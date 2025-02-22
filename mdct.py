"""
Music 422 Marina Bosi

- mdct.py -- Computes a reasonably fast MDCT/IMDCT using the FFT/IFFT

-----------------------------------------------------------------------
© 2009-2025 Marina Bosi & Richard E. Goldberg -- All rights reserved
-----------------------------------------------------------------------

"""

### ADD YOUR CODE AT THE SPECIFIED LOCATIONS ###

import numpy as np
import scipy.fft
import timeit

### Problem 1.a ###
def MDCTslow(data, a, b, isInverse=False):
    """
    Slow MDCT algorithm for window length a+b following pp. 130 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    and where the 2/N factor is included in the forward transform instead of inverse.
    a: left half-window length
    b: right half-window length
    """

    ### YOUR CODE STARTS HERE ###
    N = a + b
    n0 = (b + 1) / 2
    
    if not isInverse:
        X = np.zeros(N//2)
        for k in range(N//2):
            X[k] = (2 / N) * sum(data[n] * np.cos((2 * np.pi / N) * (n + n0) * (k + 0.5)) for n in range(N))
        return X
    else:
        x = np.zeros(N)
        for n in range(N):
            x[n] = 2 * sum(data[k] * np.cos((2 * np.pi / N) * (n + n0) * (k + 0.5)) for k in range(N//2))
        return x
    ### YOUR CODE ENDS HERE ###

### Problem 1.c ###
def MDCT(data, a, b, isInverse=False):
    """
    Fast MDCT algorithm for window length a+b following pp. 141-143 of
    Bosi & Goldberg, "Introduction to Digital Audio..." book
    and where the 2/N factor is included in forward transform instead of inverse.
    a: left half-window length
    b: right half-window length
    """

    ### YOUR CODE STARTS HERE ###
    N = a + b
    n0 = (b + 1) / 2

    if not isInverse:

        y = data * np.exp(-1j * 2 * np.pi * np.arange(N) / (2 * N))
        Y = scipy.fft.fft(y)
        X = (2 / N) * np.real(Y[:N//2] * np.exp(-1j * 2 * np.pi * n0 * (np.arange(N//2) + 0.5) / N))

        return X
    else:
        Z = np.zeros(N, dtype=complex)
        data_k = np.concatenate((data, -data[::-1]))
        Z = data_k * np.exp(1j * 2 * np.pi * np.arange(N) * n0 / N)

        z = N * scipy.fft.ifft(Z)

        x = np.real(z * np.exp(1j * 2 * np.pi * (np.arange(N) + n0) / (2 * N)))

        return x

    ### YOUR CODE ENDS HERE ###

def IMDCT(data,a,b):

    ### YOUR CODE STARTS HERE ###

    return MDCT(data, a, b, isInverse=True)
    ### YOUR CODE ENDS HERE ###

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    ### YOUR TESTING CODE STARTS HERE ###

    # Problem 1.a 
    print("Problem 1.a")
    x = [0, 1, 2, 3, 4, 4, 4, 4, 3, 1, -1, -3]
    a = 4
    b = 4
    N = a + b
    blocks = [
        [0, 0, 0, 0, 0, 1, 2, 3],
        [0, 1, 2, 3, 4, 4, 4, 4],
        [4, 4, 4, 4, 3, 1, -1, -3],
        [3, 1, -1, -3, 0, 0, 0, 0]
    ]

    for block in blocks:
        mdct_result = MDCTslow(block, a, b, isInverse=False)
        imdct_result = MDCTslow(mdct_result, a, b, isInverse=True)
        print("Input Block:", block)
        print("MDCT Result:", mdct_result)
        print("IMDCT Result:", imdct_result)
        print()

    # Problem 1.b
    print("Problem 1.b")
    x = [0, 1, 2, 3, 4, 4, 4, 4, 3, 1, -1, -3]
    a = 4
    b = 4
    N = a + b  

    output_signal = np.zeros(a)
    prior_block = [0, 0, 0, 0] 
    end_block = [0, 0, 0, 0]

    for i in range(N // 2):
        if i == (N // 2) - 1:
            new_samples = [0, 0, 0, 0]
        else:
            new_samples = x[i * b : i * b + b]

        block = prior_block + new_samples
        print("Input Block:", block)
        
        mdct = MDCTslow(block, a, b, isInverse=False)
        print("MDCT Result:", mdct)
        imdct = MDCTslow(mdct, a, b, isInverse=True) / 2
        print("IMDCT Result:", imdct)
    
        output_signal[i * b : i * b + a] += imdct[:a]
        output_signal = np.append(output_signal, imdct[-b:])
    
        prior_block = x[i * b : i * b + b]

    #output_signal = output_signal[a:-b]
    output_signal = [int(round(val)) for val in output_signal]
    print("Output Signal:", output_signal)

    # Problem 1.c 
    print("Problem 1.c")
    x = [0, 1, 2, 3, 4, 4, 4, 4, 3, 1, -1, -3]
    a = 4
    b = 4
    N = a + b
    blocks = [
        [0, 0, 0, 0, 0, 1, 2, 3],
        [0, 1, 2, 3, 4, 4, 4, 4],
        [4, 4, 4, 4, 3, 1, -1, -3],
        [3, 1, -1, -3, 0, 0, 0, 0]
    ]

    for block in blocks:
        mdct_result = MDCT(block, a, b, isInverse=False)
        imdct_result = IMDCT(mdct_result, a, b)
        print("Input Block:", block)
        print("MDCT Result:", mdct_result)
        print("IMDCT Result:", imdct_result)
        print()
    
    print("\nTiming test for block size 2048:")
    a = b = 1024 
    test_data = np.random.rand(2048) 

    slow_time = timeit.timeit(lambda: MDCTslow(test_data, a, b), number=1)
    
    fast_time = timeit.timeit(lambda: MDCT(test_data, a, b), number=1)
    
    print(f"Time for MDCTslow: {slow_time:} seconds")
    print(f"Time for MDCT: {fast_time:} seconds")
    print(f"Speedup factor: {slow_time/fast_time:.2f}x")
    
    test_freq = np.random.rand(1024)
    slow_time = timeit.timeit(lambda: MDCTslow(test_freq, a, b, isInverse=True), number=1)

    fast_time = timeit.timeit(lambda: IMDCT(test_freq, a, b), number=1)
    
    print(f"\nTime for IMDCTslow: {slow_time:} seconds")
    print(f"Time for IMDCT: {fast_time:} seconds")
    print(f"Speedup factor: {slow_time/fast_time:.2f}x")
    ### YOUR TESTING CODE ENDS HERE ###

