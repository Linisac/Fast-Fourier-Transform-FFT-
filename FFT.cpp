#include <iostream>
#include <cmath>
#include "FFT.h"

#define PI 3.141592653589793

std::complex<long double> root_of_unity(unsigned int p, int k) {
    unsigned int h = static_cast<unsigned int>(k % static_cast<int>(pow(2, p)));
    unsigned int n = static_cast<unsigned int>(pow(2, p));
    if (h == 0) return std::complex<long double>(1, 0);
    else if (4 * h == n) return std::complex<long double>(0, 1);
    else if (2 * h == n) return std::complex<long double>(-1, 0);
    else if (4 * h == 3 * n) return std::complex<long double>(0, -1);
    else return std::complex<long double>(cos(2 * h * PI / n), sin(2 * h * PI / n));
}

void f(std::vector<std::complex<long double>> &a) {
    size_t n = a.size();
    if (n > 1) {
        unsigned int power = static_cast<unsigned int>(ceil(log2(n)));
        unsigned int size = static_cast<unsigned int>(pow(2, power));
        a.resize(size);
        std::vector<std::complex<long double>> y_even(size / 2), y_odd(size / 2);
        for (unsigned int i = 0; i < size / 2; i++) y_even[i] = a[2 * i], y_odd[i] = a[2 * i + 1];
        f(y_even), f(y_odd);
        for (unsigned int i = 0; i < size / 2; i++) {
            std::complex<long double> factor = root_of_unity(power, -i);
            a[i] = y_even[i] + factor * y_odd[i];
            a[i + size / 2] = y_even[i] - factor * y_odd[i];
        }
    }
}

void FFT(std::vector<std::complex<long double>> &a) {
    size_t n = a.size();
    if (n > 1) {
        unsigned int power = static_cast<unsigned int>(ceil(log2(n)));
        unsigned int size = static_cast<unsigned int>(pow(2, power));
        a.resize(size);
        std::vector<std::complex<long double>> y_even(size / 2), y_odd(size / 2);
        for (unsigned int i = 0; i < size / 2; i++) y_even[i] = a[2 * i], y_odd[i] = a[2 * i + 1];
        FFT(y_even), FFT(y_odd);
        for (unsigned int i = 0; i < size / 2; i++) {
            std::complex<long double> factor = root_of_unity(power, i);
            a[i] = y_even[i] + factor * y_odd[i];
            a[i + size / 2] = y_even[i] - factor * y_odd[i];
        }
    }
}

void IFFT(std::vector<std::complex<long double>> &a) {
    f(a);
    unsigned int size = a.size();
    for (unsigned int i = 0; i < size; i++) a[i] /= size;
}