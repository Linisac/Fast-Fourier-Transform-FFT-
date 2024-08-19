There are two implementations of fast Fourier transform (FFT) and its inverse operation (IFFT), one in Python (FFT.py) and the other in C++11 (FFT.h and FFT.cpp).
1. Both implementations take a list/vector of complex numbers and extend its length to the closest power of 2 (if the original length is not an esact power of 2) by padding it with 0's before performing the essential operations of FFT or IFFT.
2. There is an issue of precision in the C++11 implementation.
