#ifndef FFT_LIBRARY_H
#define FFT_LIBRARY_H

#include "../../types.h"

// Complex number structure
typedef struct
{
    double real;
    double imag;
} complex_number;

int buffer_index_plus_fftLib(int buffer_next_i, int plus, int max);

int my_log2_fftLib(int N);

int reverse_fftLib(int N, int n);

void ordina_fftLib(complex_number *f1, int N);

void transform_fftLib(complex_number *f, int N);

void FFT_fftLib(complex_number *f, int N, double d);

#endif