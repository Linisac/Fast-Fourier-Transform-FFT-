import math

def FFT(a): #fast Fourier transform
    #Input:  a list 'a' of complex numbers
    #Output: the discrete Fourier transform (DFT) of 'a', where the length is the closest power of 2 >= the original length of 'a'
    n = len(a)
    if n > 1:
        a.extend([0 for _ in range(2 ** math.ceil(math.log2(n)) - n)]) #extend the 'a' so that its length is an exact power of 2
        n = len(a) #update 'n' to the new length of 'a'
        y_even = [a[i] for i in range(n) if i % 2 == 0] #even-indexed elements of 'a'
        y_odd  = [a[i] for i in range(n) if i % 2 == 1] #odd-indexed elements of 'a'
        FFT(y_even)
        FFT(y_odd)
        power = int(math.log2(n))
        for k in range(n // 2):
            factor = root_of_unity(power, k)
            a[k] = y_even[k] + factor * y_odd[k]
            a[k + n // 2] = y_even[k] - factor * y_odd[k]

def IFFT(a): #inverse fast Fourier transform
    #Input:  a list 'a' of complex numbers
    #Output: the inverse discrete Fourier-transform of 'a', where the length is the closest power of 2 >= the original length of 'a'
    def f(a): #auxiliary function that is analogous to FFT
        n = len(a)
        if n > 1:
            a.extend([0 for _ in range(2 ** math.ceil(math.log2(n)) - n)])
            n = len(a)
            y_even = [a[i] for i in range(n) if i % 2 == 0]
            y_odd  = [a[i] for i in range(n) if i % 2 == 1]
            f(y_even)
            f(y_odd)
            power = int(math.log2(n))
            for k in range(n // 2):
                factor = root_of_unity(power, -k)
                a[k] = y_even[k] + factor * y_odd[k]
                a[k + n // 2] = y_even[k] - factor * y_odd[k]
    f(a)
    n = len(a)
    for i in range(n):
        a[i] /= n

def root_of_unity(p, k):
    #Input:  integers 'p' (assumed to be > 0) and 'k'
    #Output: the kth power of the principal (2**p)th root of unity
    k = k % (2 ** p)
    if k == 0:
        return 1 + 0j
    elif 4 * k == 2 ** p:
        return 0 + 1j
    elif 2 * k == 2 ** p:
        return -1 + 0j
    elif 4 * k == 3 * (2 ** p):
        return 0 - 1j
    else:
        return math.cos(2 * k * math.pi / 2 ** p) + 1j * math.sin(2 * k * math.pi / 2 ** p)
