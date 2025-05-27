/* ----------------------------------------------------
* FFT library
* ----------------------------------------------------
* 
* This library includes the function needed to perform the FFT.
*
*/

#include "fft_library.h"
#include <math.h>
#include <stdio.h> 
#include <string.h> //I don't think it's needed 
#include <stdlib.h>
#include <stdbool.h> //I dont think it's needed

#define M_PI 3.14159265358979323846 // pi
#define WINDOW_LEN_MAX 512
#define WINDOW_LEN 128
#define N_PAD 1024

int buffer_index_plus_fftLib(int buffer_next_i, int plus, int max)
{
    return (buffer_next_i + plus) % max;
}

// Log base 2 function: Compute del log2 of N, returning the max exponent k such that 2^k <= N
int my_log2_fftLib(int N)
{
    int k = N, i = 0;
    while (k)
    {
        k >>= 1; //Shift by 1 to divide by 2
        i++;
    }
    return i - 1;
}

// Function to calculate the reverse index based on bit permutation: bit reversal operation
int reverse_fftLib(int N, int n)
{
    int log2N = my_log2_fftLib(N);
    int j, p = 0;
    for (j = 1; j <= log2N; j++)
    {
        if (n & (1 << (log2N - j)))
            p |= 1 << (j - 1);
    }
    return p;
}


/*
// Function to reorder the array based on the reverse index; 
complex_number f2[WINDOW_LEN_MAX];

void ordina_fftLib(complex_number *f1, int N)
{
    for (int i = 0; i < N; i++)
        f2[i] = f1[reverse_fftLib(N, i)];
    for (int j = 0; j < N; j++)
        f1[j] = f2[j];
}
*/

void ordina_fftLib(complex_number *f1, int N)
{
    if (N == WINDOW_LEN){
        complex_number f2[WINDOW_LEN];
        for (int i = 0; i < N; i++)
            f2[i] = f1[reverse_fftLib(N, i)];
        for (int j = 0; j < N; j++)
            f1[j] = f2[j];
    } 
    else if (N == N_PAD) {
        complex_number f2[N_PAD];
        for (int i = 0; i < N; i++)
            f2[i] = f1[reverse_fftLib(N, i)];
        for (int j = 0; j < N; j++)
            f1[j] = f2[j];
    }
    else {
        printf("Error: ordina_fftLib called with unsupported N value: %d\n", N);
    }
}

void transform_fftLib(complex_number *f, int N)
{
    ordina_fftLib(f, N); // First reorder the array
    complex_number *W;
    W = (complex_number *)malloc(N / 2 * sizeof(complex_number));
    W[1].real = cos(-2. * M_PI / N);
    W[1].imag = sin(-2. * M_PI / N);
    W[0].real = 1;
    W[0].imag = 0;
    for (int i = 2; i < N / 2; i++)
    {
        W[i].real = cos(-2. * M_PI * i / N);
        W[i].imag = sin(-2. * M_PI * i / N);
    }
    
   // It works with this:
   int step = 1;
    while (step < N) {
        int jump = step * 2;
        for (int i = 0; i < N; i += jump) {
            for (int j = 0; j < step; j++) {
                int index1 = i + j;
                int index2 = i + j + step;
                int W_index = (j * (N / jump)) % (N / 2);
                
                complex_number temp;
                temp.real = W[W_index].real * f[index2].real - W[W_index].imag * f[index2].imag;
                temp.imag = W[W_index].real * f[index2].imag + W[W_index].imag * f[index2].real;
                
                f[index2].real = f[index1].real - temp.real;
                f[index2].imag = f[index1].imag - temp.imag;
                f[index1].real += temp.real;
                f[index1].imag += temp.imag;
            }
        }
        step *= 2;
    }
    //
    free(W);
}

// FFT function
void FFT_fftLib(complex_number *f, int N, double d)
{
    transform_fftLib(f, N);
    // Scale the FFT result by multiplying each value by the step size 'd'
    for (int i = 0; i < N; i++)
    {
        f[i].real *= d;
        f[i].imag *= d;
    }
}
