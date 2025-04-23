/* ----------------------------------------------------
* AUTOCORREALATION2 ALGORITHM
* ----------------------------------------------------
* Description:
* The algorithm takes overlapping windows of 5.12s (128 samples @ 25 Hz) with step of 2.56s (64 samples @ 25Hz) and performs the autocorrelation. The autocorrelatio is
* computed using the Fourier Transform of the signal, as the IFFT of the power spectrum.
* To change the window length -> WINDOW_LEN (It must always be a power of 2). To change the overlapping -> WINDOW_STEP   
* The algorithm searches the second peak of the autocorrelation, limiting the search between the [8, 40] lags,
* corresponding to the range of HR [37.5bpm, 180bpm] or [0.625Hz, 3Hz] -> To change this parameters: NUM_AUTOCORR_LAGS, FIRST_AUTOCORR_PEAK_LAG. Note that these 
* parameters are indipendent from the window length of the signal that is considered (N of the FFT in this case).
* The algorithm finds the peak with a peak detection on the original autocorrelation signal. 
* The algorithm performs a median filter on the last 8 output -> HRM_HIST_LEN and HRM_MEDIAN_LEN
* 
* The algorithm inspired by the algorithm in "https://github.com/VirginiaSek/Algo" used for step counting from accelerometer signal.
*
* Modifications have been made in order to adapt the algorithm to the HR computation from PPG signal.
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>


#include "../../utils/bandpass_filter/bandpass_filter.h"
#include "../../types.h"


#include "autocorrelation2HeartRate.h"

#define WINDOW_LEN 128  // sliding window length, better if power of 2 (if we want to switch to FFT), 64 samples = 2.56s , 128 samples = 5.12 s, 256 samples = 10.24s
#define WINDOW_STEP 64  // step of the sliding window, 64 samples = 2.56s
#define SAMPLING_FREQ 25 // sampling frequency of the PPG signal
#define M_PI 3.14159265358979323846 // pi
//Autocorrelation lags
#define NUM_AUTOCORR_LAGS  41      // number of lags to calculate for autocorrelation. At a minimum heart rate of 37.5 bpm -> 0.625 beat/s -> 1.6 / sampling_period(0.04s) = 1.6s * 25Hz = 40 lags 
#define FIRST_AUTOCORR_PEAK_LAG 8 // corresponds to the first feasible autocorrelation lag -> at a max heart rate of 3 beats /s -> 0.333s / sampling_period(0.04s) = 8.25 -> 8 lags


// Buffers and counters
static ppg_t signal_buffer[WINDOW_LEN] = {0};
static int signal_buffer_next_i = 0;
static int HR = 0; 
static int samples_since_last_HR = 0;

//DUMP FILE: to save .csv files of FFT signals for each window
// #define DUMP_FILE 
#ifdef DUMP_FILE
static int aut2_passes = 0; // counter of how many times the autocorr has been called
#define DUMP_AUT2_FILE_NAME "aut2"
static FILE *aut2File;
#endif

//For the median filter
//Buffer per implementare media mobile o mediana sui risultati
#define HRM_HIST_LEN 8
#define HRM_MEDIAN_LEN 4
static int aut2_results[HRM_HIST_LEN] = {0};
static int aut2_results_index = 0;

//Band Pass Filter
static BPFilter bpFilter;

// Complex number structure
typedef struct
{
    double real;
    double imag;
} complex_double;

static complex_double fft_input[WINDOW_LEN];
static complex_double power_spectrum[WINDOW_LEN];


static int buffer_index_plus2(int buffer_next_i, int plus, int max)
{
    return (buffer_next_i + plus) % max;
}

// Log base 2 function: Compute del log2 of N, returning the max exponent k such that 2^k <= N
int my_log2_aut(int N)
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
int reverse_aut(int N, int n)
{
    int log2N = my_log2_aut(N);
    int j, p = 0;
    for (j = 1; j <= log2N; j++)
    {
        if (n & (1 << (log2N - j)))
            p |= 1 << (j - 1);
    }
    return p;
}

// Function to reorder the array based on the reverse index
void ordina_aut(complex_double *f1, int N)
{
    complex_double f2_1[WINDOW_LEN];
    
    for (int i = 0; i < N; i++)
        f2_1[i] = f1[reverse_aut(N, i)];
    for (int j = 0; j < N; j++)
        f1[j] = f2_1[j];
    
}

void transform_aut(complex_double *f, int N)
{
    ordina_aut(f, N); // First reorder the array
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
void FFT_aut(complex_double *f, int N, double d)
{
    transform_aut(f, N);
    // Scale the FFT result by multiplying each value by the step size 'd'
    for (int i = 0; i < N; i++)
    {
        f[i].real *= d;
        f[i].imag *= d;
    }
}


/// Initialise step counting
void autocorrelation2_heartrate_init()
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
        aut2_results[i] = 0;;
    }
    aut2_results_index = 0;

    BPFilter_init(&bpFilter);

#ifdef DUMP_FILE
    aut2_passes = 0;
#endif

}

int autocorrelation2_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz){

      
    // Applying the filter
    BPFilter_put(&bpFilter, ppg);
    ppg_t ppg_filtered = BPFilter_get(&bpFilter);


    // Add the magnitude to the circular buffer
    signal_buffer[signal_buffer_next_i] = ppg_filtered;
    signal_buffer_next_i = (signal_buffer_next_i + 1) % WINDOW_LEN;

    samples_since_last_HR++;

    // After WINDOW_STEP samples, check if it's time to perform the analysis
    if (samples_since_last_HR >= WINDOW_STEP)
    {

#ifdef DUMP_FILE
        aut2_passes++;
        char aut2FileName[100] = DUMP_AUT2_FILE_NAME;
        char idxstr[5];
        sprintf(idxstr, "%d", aut2_passes);
        strcat(aut2FileName, idxstr);
        strcat(aut2FileName, ".csv");
        aut2File = fopen(aut2FileName, "w+");
#endif

        
        samples_since_last_HR = 0;

        // Prepare the data for FFT
        for (int i = 0; i < WINDOW_LEN; i++)
        {
            if (i < WINDOW_LEN){
                    int buffer_i = buffer_index_plus2(signal_buffer_next_i, i, WINDOW_LEN);
                    fft_input[i].real = (double)signal_buffer[buffer_i];
                    fft_input[i].imag = 0.0;
            } 
            
        }

        // Perform the FFT
        FFT_aut(fft_input, WINDOW_LEN, 1.0);



        for (int i = 0; i < WINDOW_LEN; i++) //Computed on all the frequencies
        {
            power_spectrum[i].real = fft_input[i].real * fft_input[i].real + fft_input[i].imag * fft_input[i].imag;
            power_spectrum[i].imag = 0.0;
        }

        //Implemented the IFFT in this way-> Conjugate the values, perform the FFT, conjugate again and divide by N (equivalent to do IFFT)
        // 1. Conjugate the values
        for (int i = 0; i < WINDOW_LEN; i++) {
            power_spectrum[i].imag *= -1.0; //This passage is complitely useless because the imaginary part is zero, but we still keep it!
        }
        //2. FFT 
        FFT_aut(power_spectrum, WINDOW_LEN, 1.0);
        //3. Conjugate again and divide by N
        for (int i = 0; i < WINDOW_LEN; i++) {
            power_spectrum[i].imag *= -1.0;
            power_spectrum[i].real /= WINDOW_LEN;
            power_spectrum[i].imag /= WINDOW_LEN;
        }
        
        //Write on the file
#ifdef DUMP_FILE
        for (int i=0;i<WINDOW_LEN;i++){
            if (aut2File){
                    fprintf(aut2File, "%d,%f\n", i,power_spectrum[i].real);
                }
        }

#endif


#ifdef DUMP_FILE
        if (aut2File)
        {
            fflush(aut2File);
            fclose(aut2File);
        }
#endif

        //In power spectrum there is the autocorrelation -> find the lag (index) of the peak
        uint8_t peak_ind = 0;
        float interval = 0;
        uint8_t i;
        for (i = FIRST_AUTOCORR_PEAK_LAG; i < NUM_AUTOCORR_LAGS; i++)
        {

            if ((power_spectrum[i].real > power_spectrum[i - 1].real) && (power_spectrum[i].real > power_spectrum[i + 1].real)){
                peak_ind = i;
                break; 
            }
            
        }
        if (peak_ind != 0) {
            interval = ((float)(peak_ind) / (float)(SAMPLING_FREQ)); //This is the period T in seconds of the periodic signal
        } else {
            interval = 0;
        }
     
        //If the peak has not been found (interval = 0), just keep the previous value of HR
        if (interval!=0){
            HR = (int)(600/interval); //bpm*10
        } 

        //For the median filter
        //Saving the results in a buffer to implement a moving average or median
        aut2_results[aut2_results_index] = HR; 
        aut2_results_index++;
        if (aut2_results_index >= HRM_HIST_LEN)
        {
            aut2_results_index = 0;
        }
        //Median
        //First: Sorting the buffer
        bool busy;
        do {
            busy = false;
            for (int i=0;i<HRM_HIST_LEN-1;i++) {
                if (aut2_results[i] > aut2_results[i+1]) {
                    int te = aut2_results[i];
                    aut2_results[i] = aut2_results[i+1];
                    aut2_results[i+1] = te;
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
            if (aut2_results[i]==0) continue;
            sumBPM += aut2_results[i];
            n++;
        }
        if (n) {
            HR = (int)(sumBPM/n);
        } 
    }

    // Return the HR*10
    return (int)(HR);

}