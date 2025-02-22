"""
Music 422 - Marina Bosi

quantize.py -- routines to quantize and dequantize floating point values
between -1.0 and 1.0 ("signed fractions")

-----------------------------------------------------------------------
© 2009-2025 Marina Bosi & Richard E. Goldberg -- All rights reserved
-----------------------------------------------------------------------
"""

### ADD YOUR CODE AT THE SPECIFIED LOCATIONS ###

import numpy as np

### Problem 1.a.i ###
def QuantizeUniform(aNum,nBits):
    """
    Uniformly quantize signed fraction aNum with nBits
    """
    #Notes:
    #The overload level of the quantizer should be 1.0

    ### YOUR CODE STARTS HERE ###
    aNum = np.clip(aNum, -1.0, 1.0)
    max_code = int(1 << (nBits - 1)) - 1

    code = round(aNum * max_code)
    return np.clip(code, -max_code, max_code) 
    # s = 0 if aNum >= 0 else 1
    # if abs(aNum) >= 1:
    #     code = 2 ** (nBits - 1) - 1
    # else:
    #     code = int(((2 ** nBits - 1) * abs(aNum) + 1) / 2)

    # aQuantizedNum = ((-1) ** s) * abs(code)
    ### YOUR CODE ENDS HERE ###
 
    return aQuantizedNum

### Problem 1.a.i ###
def DequantizeUniform(aQuantizedNum,nBits):
    """
    Uniformly dequantizes nBits-long number aQuantizedNum into a signed fraction
    """

    ### YOUR CODE STARTS HERE ###
    max_code = int(1 << (nBits - 1)) - 1
    return np.clip(aQuantizedNum / max_code, -1.0, 1.0)
    # sign = 1 if aQuantizedNum >= 0 else -1
    # number = 2.0 * abs(aQuantizedNum) / (2 ** nBits - 1)
    # aNum = sign * abs(number)
    ### YOUR CODE ENDS HERE ###

    return aNum

### Problem 1.a.ii ###
def vQuantizeUniform(aNumVec, nBits):
    """
    Uniformly quantize vector aNumberVec of signed fractions with nBits
    """

    #Notes:
    #Make sure to vectorize properly your function as specified in the homework instructions

    ### YOUR CODE STARTS HERE ###
    aNumVec = np.clip(aNumVec, -1.0, 1.0)
    max_code = int(1 << int(nBits - 1)) - 1
    codes = np.round(aNumVec * max_code)
    return np.clip(codes, -max_code, max_code)
    # aNumVec = np.array(aNumVec)
    # s = np.where(aNumVec >= 0, 0, 1)
    # codes = np.where(np.abs(aNumVec) >= 1, 2**(nBits - 1) - 1, np.floor(((2 ** nBits - 1) * np.abs(aNumVec) + 1) / 2))
    # aQuantizedNumVec = np.where(s == 0, codes, -codes)
    ### YOUR CODE ENDS HERE ###

    return aQuantizedNumVec

### Problem 1.a.ii ###
def vDequantizeUniform(aQuantizedNumVec, nBits):
    """
    Uniformly dequantizes vector of nBits-long numbers aQuantizedNumVec into vector of  signed fractions
    """

    ### YOUR CODE STARTS HERE ###
    max_code = int(1 << int(nBits - 1)) - 1
    return np.clip(aQuantizedNumVec / max_code, -1.0, 1.0)
    # aQuantizedNumVec = np.array(aQuantizedNumVec)
    # s = np.sign(aQuantizedNumVec)
    # numbers = 2.0 * np.abs(aQuantizedNumVec) / (2 ** nBits - 1)
    # aNumVec = s * numbers
    ### YOUR CODE ENDS HERE ###

    return aNumVec

### Problem 1.b ###
def ScaleFactor(aNum, nScaleBits=3, nMantBits=5):
    """
    Return the floating-point scale factor for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """
    #Notes:
    #The scale factor should be the number of leading zeros

    ### YOUR CODE STARTS HERE ###
    if aNum == 0:
        return (1 << nScaleBits) - 1
    
    R = (1 << nScaleBits) - 1 + nMantBits
    quantized = abs(QuantizeUniform(aNum, R))
    binary_rep = bin(quantized)[2:].zfill(R)
    leading_zeros = len(binary_rep) - len(binary_rep.lstrip('0'))
    
    return min(leading_zeros, (1 << nScaleBits) - 1)
    # R = 2**nScaleBits - 1 + nMantBits
    # quantized = QuantizeUniform(aNum, R)
    # binary = bin(abs(int(quantized)))[2:].zfill(R)
    # leading_zeros = len(binary) - len(binary.lstrip('0'))
    # scale = min(leading_zeros, 2**nScaleBits - 1)
    ### YOUR CODE ENDS HERE ###

    return scale

### Problem 1.b ###
def MantissaFP(aNum, scale, nScaleBits=3, nMantBits=5):
    """
    Return the floating-point mantissa for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """

    ### YOUR CODE STARTS HERE ###
    R = (1 << nScaleBits) - 1 + nMantBits
    quantized = abs(QuantizeUniform(aNum, R))
    sign_bit = 1 if aNum < 0 else 0
    binary_rep = bin(quantized)[2:].zfill(R)

    if scale < (1 << nScaleBits) - 1:
        mantissa_bits = binary_rep[scale+1:scale+nMantBits]
    else:
        mantissa_bits = binary_rep[scale:scale+nMantBits-1]

    mantissa = int(str(sign_bit) + mantissa_bits, 2)
    return np.clip(mantissa, 0, (1 << nMantBits) - 1)
    # R = 2**nScaleBits - 1 + nMantBits
    # quantized = QuantizeUniform(abs(aNum), R)
    # s = 0 if aNum >= 0 else 1

    # binary = bin(int(quantized))[2:].zfill(R)
    
    # if scale == 2**nScaleBits - 1:
    #     mantissa_bits = binary[scale:scale+nMantBits-1]
    # else:
    #     mantissa_bits = binary[scale+1:scale+nMantBits]

    # mantissa_bits = str(s) + mantissa_bits
    # mantissa = int(mantissa_bits, 2)
    #print(f"Input: {aNum}, Scale: {scale}, Mantissa Bits: {mantissa_bits}")

    ### YOUR CODE ENDS HERE ###

    return mantissa

### Problem 1.b ###
def DequantizeFP(scale, mantissa, nScaleBits=3, nMantBits=5):
    """
    Returns a  signed fraction for floating-point scale and mantissa given specified scale and mantissa bits
    """

    ### YOUR CODE STARTS HERE ###
    R = (1 << nScaleBits) - 1 + nMantBits
    mantissa = np.clip(mantissa, 0, (1 << nMantBits) - 1)

    sign = (mantissa >> (nMantBits - 1)) & 1
    mantissa_bits = bin(mantissa & ((1 << (nMantBits - 1)) - 1))[2:].zfill(nMantBits - 1)

    if scale < (1 << nScaleBits) - 1:
        code_bits = '0' * scale + '1' + mantissa_bits
    else:
        code_bits = '0' * scale + mantissa_bits

    code_bits = code_bits.ljust(R, '0')[:R]
    code = int(code_bits, 2)
    
    return -DequantizeUniform(code, R) if sign else DequantizeUniform(code, R)
    # R = 2**nScaleBits - 1 + nMantBits

    # s = (mantissa >> (nMantBits - 1)) & 1

    # mantissa_bits = bin(mantissa & ((1 << (nMantBits - 1)) - 1))[2:].zfill(nMantBits - 1)

    # if scale == 2**nScaleBits - 1:
    #     code_bits = '0' * scale + mantissa_bits
    # else:
    #     code_bits = '0' * scale + '1' + mantissa_bits

    # code_bits = code_bits = code_bits.ljust(R, '0')[:R]

    # code = int(code_bits, 2)

    # if s == 1:
    #     code = -code 

    # aNum = DequantizeUniform(code, R)
    #print(f"Code Bits: {code_bits}, Code: {code}, Reconstructed: {aNum}")
    ### YOUR CODE ENDS HERE ###

    return aNum

### Problem 1.c.i ###
def Mantissa(aNum, scale, nScaleBits=3, nMantBits=5):
    """
    Return the block floating-point mantissa for a  signed fraction aNum given nScaleBits scale bits and nMantBits mantissa bits
    """

    ### YOUR CODE STARTS HERE ###
    R = (1 << nScaleBits) - 1 + nMantBits
    quantized = abs(QuantizeUniform(aNum, R))
    sign_bit = 1 if aNum < 0 else 0
    binary_rep = bin(quantized)[2:].zfill(R)

    mantissa_bits = binary_rep[scale:scale + nMantBits - 1]
    mantissa = int(str(sign_bit) + mantissa_bits, 2)

    return np.clip(mantissa, 0, int(1 << nMantBits) - 1)
    # R = 2**nScaleBits - 1 + nMantBits 
    # quantized = QuantizeUniform(abs(aNum), R) 
    # s = 0 if aNum >= 0 else 1  

    # binary = bin(int(quantized))[2:].zfill(R)

    # mantissa_bits = binary[scale:scale + nMantBits - 1]
    # mantissa_bits = str(s) + mantissa_bits
    # mantissa = int(mantissa_bits, 2)
    ### YOUR CODE ENDS HERE ###

    return mantissa

### Problem 1.c.i ###
def Dequantize(scale, mantissa, nScaleBits=3, nMantBits=5):
    """
    Returns a  signed fraction for block floating-point scale and mantissa given specified scale and mantissa bits
    """

    ### YOUR CODE STARTS HERE ###
    R = (1 << nScaleBits) - 1 + nMantBits
    mantissa = np.clip(mantissa, 0, int(1 << nMantBits) - 1)

    sign = (mantissa >> int(nMantBits - 1)) & 1
    mantissa_bits = bin(mantissa & ((1 << int(nMantBits - 1)) - 1))[2:].zfill(nMantBits - 1)

    code_bits = '0' * scale + mantissa_bits
    code_bits = code_bits.ljust(R, '0')[:R]
    code = int(code_bits, 2)
    
    return -DequantizeUniform(code, R) if sign else DequantizeUniform(code, R)
    # R = 2**nScaleBits - 1 + nMantBits 
    # s = (mantissa >> (nMantBits - 1)) & 1
    # mantissa_bits = bin(mantissa & ((1 << (nMantBits - 1)) - 1))[2:].zfill(nMantBits - 1)

    # if int(mantissa_bits, 2) == 0:
    #     code_bits = '0' * scale + mantissa_bits
    # else:
    #     code_bits = '0' * scale + mantissa_bits

    # code_bits = code_bits.ljust(R, '0')[:R]
    # code = int(code_bits, 2)
    # if s == 1:
    #     code = -code 

    # aNum = DequantizeUniform(code, R)
    ### YOUR CODE ENDS HERE ###

    return aNum

### Problem 1.c.ii ###
def vMantissa(aNumVec, scale, nScaleBits=3, nMantBits=5):
    """
    Return a vector of block floating-point mantissas for a vector of  signed fractions aNum given nScaleBits scale bits and nMantBits mantissa bits
    """

    ### YOUR CODE STARTS HERE ###
    R = 2**nScaleBits - 1 + nMantBits
    quantizedVec = np.abs(vQuantizeUniform(aNumVec, R))
    s = np.where(aNumVec >= 0, 0, 1)
    binaryVec = np.array([bin(int(q))[2:].zfill(R) for q in quantizedVec])
    mantissa_bits = np.array([b[scale:scale + nMantBits - 1] for b in binaryVec])
    mantissaVec = s * (1 << (nMantBits - 1)) + np.array([int(m.strip(), 2) for m in mantissa_bits], dtype=int)
    ### YOUR CODE ENDS HERE ###

    return mantissaVec

### Problem 1.c.ii ###
def vDequantize(scale, mantissaVec, nScaleBits=3, nMantBits=5):
    """
    Returns a vector of  signed fractions for block floating-point scale and vector of block floating-point mantissas given specified scale and mantissa bits
    """

    ### YOUR CODE STARTS HERE ###
    R = 2**nScaleBits - 1 + nMantBits

    sign_bits = (mantissaVec >> (nMantBits - 1)) & 1
    signs = np.where(sign_bits == 1, -1, 1)

    mantissa_bits = mantissaVec & (int(1 << int(nMantBits - 1)) - 1)

    binaryVec = np.array([bin(int(m))[2:].zfill(nMantBits - 1) for m in mantissa_bits])
    code_bits = np.array(['0' * scale + b.ljust(R - scale, '0') for b in binaryVec])
    codes = np.array([int(c, 2) for c in code_bits])

    aNumVec = signs * vDequantizeUniform(codes, R)
    ### YOUR CODE ENDS HERE ###

    return aNumVec

#-----------------------------------------------------------------------------

#Testing code
if __name__ == "__main__":

    ### YOUR TESTING CODE STARTS HERE ###

    inputs = [-0.99, -0.38, -0.10, -0.01, -0.001, 0.0, 0.05, 0.28, 0.65, 0.97, 1.0]

    # 8 bit midtread test
    outputs = []
    print("Testing 8 bit uniform midtread quantization...")
    for input in inputs:
        outputs.append(DequantizeUniform(QuantizeUniform(input, 8), 8))
    voutputs = vDequantizeUniform(vQuantizeUniform(inputs, 8), 8)
    for i in range(len(inputs)):
        print(f"Input: {inputs[i]}, Output: {outputs[i]:.4f}, vOutput: {voutputs[i]:.4f}")

    # 12 bit midtread test
    outputs = []
    print("Testing 12 bit uniform midtread quantization...")
    for input in inputs:
        outputs.append(DequantizeUniform(QuantizeUniform(input, 12), 12))
    voutputs = vDequantizeUniform(vQuantizeUniform(inputs, 12), 12)
    for i in range(len(inputs)):
        print(f"Input: {inputs[i]}, Output: {outputs[i]:.4f}, vOutput: {voutputs[i]:.4f}")

    # Floating-point quantization test
    outputs = []
    nScaleBits = 3
    nMantBits = 5
    print("Testing floating-point quantization...")
    for input in inputs:
        scale = ScaleFactor(input, nScaleBits, nMantBits)
        mantissa = MantissaFP(input, scale, nScaleBits, nMantBits)
        outputs.append(DequantizeFP(scale, mantissa, nScaleBits, nMantBits))
    
    for i in range(len(inputs)):
        print(f"Input: {inputs[i]}, Output: {outputs[i]:.4f}")

    # Block floating-point quantization test
    outputs = []
    nScaleBits = 3
    nMantBits = 5
    print("Testing block floating-point quantization...")
    for input in inputs:
        scale = ScaleFactor(input, nScaleBits, nMantBits)
        mantissa = Mantissa(input, scale, nScaleBits, nMantBits)
        outputs.append(Dequantize(scale, mantissa, nScaleBits, nMantBits))
    
    for i in range(len(inputs)):
        print(f"Input: {inputs[i]}, Output: {outputs[i]:.4f}")


    print("Testing vectorized block floating-point quantization...")
    aNumVec = np.array(inputs)
    nScaleBits = 3
    nMantBits = 5
    scale = ScaleFactor(np.max(np.abs(aNumVec)), nScaleBits, nMantBits)

    mantissaVec = vMantissa(aNumVec, scale, nScaleBits, nMantBits)
    aNumVec_dequantized = vDequantize(scale, mantissaVec, nScaleBits, nMantBits)

    for i in range(len(aNumVec)):
        print(f"Input: {aNumVec[i]}, Output: {aNumVec_dequantized[i]:.4f}")
    ### YOUR TESTING CODE ENDS HERE ###



