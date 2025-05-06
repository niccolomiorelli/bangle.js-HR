/* ----------------------------------------------------
* FFT ALGORITHM
* ----------------------------------------------------
* Description:
* The algorithm takes overlapping windows of 5.12s (128 samples @ 25 Hz) with step of 2.56s (64 samples @ 25Hz) and applies the Fast Fourier Transform.
* To change the window length -> WINDOW_LEN. To change the overlapping -> WINDOW_STEP   
* Note: if the WINDOW_len is changed, the MIN_FREQ_FFT_I and MAX_FREQ_FFT_I must be changed too.
* The algorithm searches the peak of the Fourier transform, limiting the search between the indexes [3, 22],
* corresponding to the range of HR [35bpm, 246 bpm] or [0.58Hz, 4.29Hz] -> To change this parameters: MIN_FREQ_FFT_I, MAX_FREQ_FFT_I. Formula: freq = i*f_sample/N
* The algorithm performs a median filter on the last 8 output -> HRM_HIST_LEN and HRM_MEDIAN_LEN
* 
* The algorithm inspired by the algorithm in "https://github.com/VirginiaSek/Algo" used for step counting from accelerometer signal.
*
* Modifications have been made in order to adapt the algorithm to the HR computation from PPG signal.
*
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>


#include "../../utils/bandpass_filter/bandpass_filter.h"
#include "../../types.h"

#include "fftHeartRate.h"

#define WINDOW_LEN 128  // sliding window length, better if power of 2 (if we want to switch to FFT), 64 samples = 2.56s , 128 samples = 5.12 s, 256 samples = 10.24s
#define WINDOW_STEP 64 // step of the sliding window, 64 samples = 2.56s
#define SAMPLING_FREQ 25 // sampling frequency of the PPG signal
#define MIN_FREQ_FFT_I 4 // index of the FFT corresponding to the minimum heart rate -> 3: corresponds to 35 bpm (if N=256  -> 5 corresponds to 30bpm)
#define MAX_FREQ_FFT_I 22 // index of the FFT corresponding to the maximum heart rate -> 22: correspnds to 246 bpm (if N=256  -> 43 corresponds to 258bpm)
#define M_PI 3.14159265358979323846 // pi

// Buffers and counters
static ppg_t signal_buffer[WINDOW_LEN] = {0};
static int signal_buffer_next_i = 0;
static int HR = 0; 
static int samples_since_last_HR = 0;

//DUMP FILE: to save .csv files of FFT signals for each window
//#define DUMP_FILE 
#ifdef DUMP_FILE
static int fft_passes = 0; // counter of how many times the autocorr has been called
#define DUMP_FFT_FILE_NAME "fft"
static FILE *fftFile;
#endif

//For the median filter
//Buffer per implementare media mobile o mediana sui risultati
#define HRM_HIST_LEN 8
#define HRM_MEDIAN_LEN 4
static int fft_results[HRM_HIST_LEN] = {0};
static int fft_results_index = 0;

//Band Pass Filter
static BPFilter bpFilter;

// Complex number structure
typedef struct
{
    double real;
    double imag;
} complex_double;

static complex_double fft_input[WINDOW_LEN];

static int buffer_index_plus(int buffer_next_i, int plus, int max)
{
    return (buffer_next_i + plus) % max;
}

// Log base 2 function: Compute del log2 of N, returning the max exponent k such that 2^k <= N
int my_log2(int N)
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
int reverse(int N, int n)
{
    int log2N = my_log2(N);
    int j, p = 0;
    for (j = 1; j <= log2N; j++)
    {
        if (n & (1 << (log2N - j)))
            p |= 1 << (j - 1);
    }
    return p;
}

// Function to reorder the array based on the reverse index
void ordina(complex_double *f1, int N)
{
    complex_double f2[WINDOW_LEN];
    for (int i = 0; i < N; i++)
        f2[i] = f1[reverse(N, i)];
    for (int j = 0; j < N; j++)
        f1[j] = f2[j];
}

void transform(complex_double *f, int N)
{
    ordina(f, N); // First reorder the array
    complex_double *W;
    W = (complex_double *)malloc(N / 2 * sizeof(complex_double));
    W[1].real = cos(-2. * M_PI / N);
    W[1].imag = sin(-2. * M_PI / N);
    W[0].real = 1;
    W[0].imag = 0;
    for (int i = 2; i < N / 2; i++)
    {
        W[i].real = cos(-2. * M_PI * i / N);
        W[i].imag = sin(-2. * M_PI * i / N);
    }
    /* //DA COMMENTARE QUESTA PARTE
    int n = 1;
    int a = N / 2;
    for (int j = 0; j < my_log2(N); j++)
    {
        for (int i = 0; i < N; i++)
        {
            if (!(i & n))
            {
                complex_double temp = f[i];
                complex_double Temp = W[(i * a) % (n * a)];
                Temp.real *= f[i + n].real - f[i].real;
                Temp.imag *= f[i + n].imag - f[i].imag;
                f[i].real += Temp.real;
                f[i].imag += Temp.imag;
                f[i + n].real = temp.real - Temp.real;
                f[i + n].imag = temp.imag - Temp.imag;
            }
        }
        n *= 2;
        a = a / 2;
    }
    */ //FINE COMMENTO DI QUESTA PARTE
   // It works with this:
   int step = 1;
    while (step < N) {
        int jump = step * 2;
        for (int i = 0; i < N; i += jump) {
            for (int j = 0; j < step; j++) {
                int index1 = i + j;
                int index2 = i + j + step;
                int W_index = (j * (N / jump)) % (N / 2);
                
                complex_double temp;
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
void FFT(complex_double *f, int N, double d)
{
    transform(f, N);
    // Scale the FFT result by multiplying each value by the step size 'd'
    for (int i = 0; i < N; i++)
    {
        f[i].real *= d;
        f[i].imag *= d;
    }
}


/// Initialise step counting
void fft_heartrate_init()
{
    HR = 0;
    samples_since_last_HR = 0;
    // Initialize the signal buffer to zeros
    for (int i = 0; i < WINDOW_LEN; i++)
    {
        signal_buffer[i] = 0;
    }
    signal_buffer_next_i = 0;

    for (int i = 0; i < HRM_HIST_LEN; i++)
    {
        fft_results[i] = 0;;
    }
    fft_results_index = 0;

    BPFilter_init(&bpFilter);

#ifdef DUMP_FILE
    fft_passes = 0;
#endif

}

int fft_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz){

      
    // Applying the filter
    BPFilter_put(&bpFilter, ppg);
    ppg_t ppg_filtered = BPFilter_get(&bpFilter);


    // Add the magnitude to the circular buffer
    signal_buffer[signal_buffer_next_i] = ppg_filtered;
    signal_buffer_next_i = (signal_buffer_next_i + 1) % WINDOW_LEN;

    samples_since_last_HR++;

    // After WINDOW_STEP samples, check if it's time to perform the analysis
    if (samples_since_last_HR == WINDOW_STEP)
    {

#ifdef DUMP_FILE
        fft_passes++;
        char fftFileName[100] = DUMP_FFT_FILE_NAME;
        char idxstr[5];
        sprintf(idxstr, "%d", fft_passes);
        strcat(fftFileName, idxstr);
        strcat(fftFileName, ".csv");
        fftFile = fopen(fftFileName, "w+");
#endif

        
        samples_since_last_HR = 0;

        // Prepare the data for FFT
        for (int i = 0; i < WINDOW_LEN; i++)
        {
            int buffer_i = buffer_index_plus(signal_buffer_next_i, i, WINDOW_LEN);
            fft_input[i].real = (double)signal_buffer[buffer_i];
            fft_input[i].imag = 0.0;
        }

        // Perform the FFT
        FFT(fft_input, WINDOW_LEN, 1.0);

        // Find the dominant frequency
        double max_fft_magnitude = 0.0;
        int dominant_freq_index = 0;

        double fft_magnitude[MAX_FREQ_FFT_I] = {0.0}; //The first ones - until MIN_FREQ_FFT_I - are initialized to zero and they keep being zero

        for (int i = MIN_FREQ_FFT_I; i < MAX_FREQ_FFT_I; i++)
        {
            //NOTE: With the following part commented: the magintude signal is saved entirely so it can be saved in a file (DUMP FILE)
            //Otherwise, just uncomment the following, for a faster computation, but the magnitude file is not saved, only the peak is searched

            //double fft_magnitude = sqrt(fft_input[i].real * fft_input[i].real + fft_input[i].imag * fft_input[i].imag);
            /*
            if (fft_magnitude > max_fft_magnitude)
            {
                max_fft_magnitude = fft_magnitude;
                dominant_freq_index = i;
            }
            */

            fft_magnitude[i] = sqrt(fft_input[i].real * fft_input[i].real + fft_input[i].imag * fft_input[i].imag);

            //Write on the file
#ifdef DUMP_FILE
            if (fftFile)
            {
                fprintf(fftFile, "%d, %f\n", i, fft_magnitude[i]);
            }
#endif
        }
        for (int i = MIN_FREQ_FFT_I; i < MAX_FREQ_FFT_I; i++)
        {
            if (fft_magnitude[i] > max_fft_magnitude)
            {
                max_fft_magnitude = fft_magnitude[i];
                dominant_freq_index = i;
            }
        }

        // Calculate the dominant frequency in Hz
        double dominant_freq = (double)dominant_freq_index * SAMPLING_FREQ / WINDOW_LEN;

        // Calculate the number of steps based on the dominant frequency
        HR  = dominant_freq * 60;

#ifdef DUMP_FILE
        if (fftFile)
        {
            fflush(fftFile);
            fclose(fftFile);
        }
#endif
        
        //For the median filter
        //Saving the results in a buffer to implement a moving average or median
        fft_results[fft_results_index] = HR; 
        fft_results_index++;
        if (fft_results_index >= HRM_HIST_LEN)
        {
            fft_results_index = 0;
        }
        //Median
        //First: Sorting the buffer
        bool busy;
        do {
            busy = false;
            for (int i=0;i<HRM_HIST_LEN-1;i++) {
                if (fft_results[i] > fft_results[i+1]) {
                    int te = fft_results[i];
                    fft_results[i] = fft_results[i+1];
                    fft_results[i+1] = te;
                    busy = true;
                }
            }
        } while (busy);
        // calculate HRM from middle values
        int16_t min = (int16_t)((HRM_HIST_LEN - HRM_MEDIAN_LEN)/2); // 2
        int16_t max = (int16_t)((HRM_HIST_LEN + HRM_MEDIAN_LEN)/2); // 6
        int16_t n = 0;
        int sumBPM = 0;
        for (int16_t i=min;i<max;i++) {
            if (fft_results[i]==0) continue;
            sumBPM += fft_results[i];
            n++;
        }
        if (n) {
            HR = (int)(sumBPM/n);
        }
        
    }

    // Return the HR*10
    return (int)(HR*10);

}