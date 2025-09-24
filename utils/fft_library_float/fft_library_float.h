#ifndef FFT_LIBRARY_FLOAT_H
#define FFT_LIBRARY_FLOAT_H

#include <stdint.h>

#ifndef WINDOW_LEN
#define WINDOW_LEN 128
#endif

#ifndef N_PAD
#define N_PAD 1024
#endif

typedef struct {
    float real;
    float imag;
} complex_number_float;

// util
int my_log2_fftLib_float(int N);
int reverse_fftLib_float(int N, int n);

// core
void FFT_fftLib_float(complex_number_float *f, float d);

// helpers
void apply_hann_window_fftLib_float(const float *input, float *windowed_output, int len);
void zero_pad_fftLib_float(const complex_number_float *in, complex_number_float *out, int N);

#endif
