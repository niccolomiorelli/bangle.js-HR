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
#include "../../utils/rolling_stats/rolling_stats.h"
#include "../../types.h"
#include "../../utils/fft_library/fft_library.h"

#include "fft2HeartRate.h"

#define WINDOW_LEN 128  // sliding window length, better if power of 2 (if we want to switch to FFT), 64 samples = 2.56s , 128 samples = 5.12 s, 256 samples = 10.24s
#define WINDOW_STEP 64 // step of the sliding window, 64 samples = 2.56s
#define SAMPLING_FREQ 25 // sampling frequency of the PPG signal
//In case of no padding:
//#define MIN_FREQ_FFT_I 4 // index of the FFT corresponding to the minimum heart rate -> 4: corresponds to 46 bpm (if N=256  -> 5 corresponds to 30bpm)
//#define MAX_FREQ_FFT_I 18 // index of the FFT corresponding to the maximum heart rate -> 18: correspnds to 210 bpm (if N=256  -> 43 corresponds to 258bpm)
//In case of padding:
#define MIN_FREQ_FFT_I 32 // 32-> with N_PAD = 1024 -> 46 bpm
#define MAX_FREQ_FFT_I 144 // 144-> with N_PAD = 1024 -> 210 bpm


#define M_PI 3.14159265358979323846 // pi

#define N_PAD 1024 //Length of the window after padding


// Buffers and counters
static float signal_buffer[WINDOW_LEN] = {0};
static int signal_buffer_next_i = 0;
static int HR = 0; 
static int samples_since_last_HR = 0;

//DUMP FILE: to save .csv files of FFT signals for each window
// #define DUMP_FILE 
#ifdef DUMP_FILE
static int fft_passes = 0; // counter of how many times the autocorr has been called
#define DUMP_FFT2_FILE_NAME "fft2"
static FILE *fft2File;
#define DUMP_ROLL_STATS_FFT2_FILE_NAME "roll_stats_fft2.csv"
static FILE *roll_stats_fft2File;
#endif

//For the median filter
//Buffer per implementare media mobile o mediana sui risultati
#define HRM_HIST_LEN 8
#define HRM_MEDIAN_LEN 4
static int fft_results[HRM_HIST_LEN] = {0};
static int fft_results_index = 0;

//Band Pass Filter
static BPFilter bpFilter;
//Standardization
static Stats stats_ppg;

static double windowed_signal_in[WINDOW_LEN]; //For the windowing
static double windowed_signal_out[WINDOW_LEN]; 
static complex_number fft_input[WINDOW_LEN];
static complex_number fft_input_padded[N_PAD];

// For the linear interpolation
typedef struct {
    float ppg;
    float accx;
    float accy;
    float accz;
    uint32_t timestamp_ms;             
} Sample;

static ppg_t ppg_last;
static accel_t accx_last;
static accel_t accy_last;
static accel_t accz_last;



// Crea la finestra di Hann
void apply_hann_window(double *input, double *windowed_output, int len) {
    for (int n = 0; n < len; n++) {
        double hann = 0.5 * (1.0 - cos(2.0 * M_PI * n / (len - 1)));
        windowed_output[n] = input[n] * hann;
    }
}


//Padding function: it adds zeros to the input, creating another complex_double vector, so I have the two input ready to be tested
void zero_pad(complex_number *in, complex_number *out, int N, int N_pad) {
    for (int i = 0; i < N_pad; i++) {
        if (i < N) {
            out[i] = in[i];
        } else {
            out[i].real = 0.0;
            out[i].imag = 0.0;
        }
    }
}


/// Initialise step counting
void fft2_heartrate_init()
{
    // Linear interpolation
    ppg_last=0;
    accx_last=0;
    accy_last=0;
    accz_last=0;    

    HR = 0;
    samples_since_last_HR = 0;
    // Initialize the signal buffer to zeros
    for (int i = 0; i < WINDOW_LEN; i++)
    {
        signal_buffer[i] = 0.0;
        windowed_signal_in[i] = 0.0;
        windowed_signal_out[i] = 0.0;
    }
    signal_buffer_next_i = 0;

    for (int i = 0; i < HRM_HIST_LEN; i++)
    {
        fft_results[i] = 0;;
    }
    fft_results_index = 0;

    BPFilter_init(&bpFilter);

    rolling_stats_reset(&stats_ppg);

#ifdef DUMP_FILE
    fft_passes = 0;
    roll_stats_fft2File = fopen(DUMP_ROLL_STATS_FFT2_FILE_NAME, "w+");
#endif


}

int main_algorithm_fft2(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz){

      
    // Applying the filter
    BPFilter_put(&bpFilter, ppg);
    ppg_t ppg_filtered = BPFilter_get(&bpFilter);

    //Standardization (optional: in case just comment it and insert ppg_filtered in the signal_buffer, instead of ppg_standardized)
    rolling_stats_addValue((float)ppg_filtered, &stats_ppg);
    float mean_ppg = rolling_stats_get_mean(&stats_ppg);
    float var_ppg = rolling_stats_get_variance(&stats_ppg);
    float std_ppg = rolling_stats_get_standard_deviation(&stats_ppg);
    if (std_ppg == 0.0) {
        std_ppg = 1.0; // Avoid division by zero
    }
    float ppg_standardized = (float)(ppg_filtered - mean_ppg) / std_ppg;

#ifdef DUMP_FILE
        if (roll_stats_fft2File)
        {
            if (!fprintf(roll_stats_fft2File, "%d, %f, %f, %f, %f\n", ppg_filtered, ppg_standardized, mean_ppg, var_ppg, std_ppg ))
                puts("error writing file");
            fflush(roll_stats_fft2File);
        }
#endif

    // Add the magnitude to the circular buffer
    signal_buffer[signal_buffer_next_i] = ppg_standardized;
    signal_buffer_next_i = (signal_buffer_next_i + 1) % WINDOW_LEN;

    samples_since_last_HR++;

    // After WINDOW_STEP samples, check if it's time to perform the analysis
    if (samples_since_last_HR == WINDOW_STEP)
    {

#ifdef DUMP_FILE
        fft_passes++;
        char fft2FileName[100] = DUMP_FFT2_FILE_NAME;
        char idxstr[5];
        sprintf(idxstr, "%d", fft_passes);
        strcat(fft2FileName, idxstr);
        strcat(fft2FileName, ".csv");
        fft2File = fopen(fft2FileName, "w+");
#endif

        
        samples_since_last_HR = 0;

        // Prepare the data for FFT
        for (int i = 0; i < WINDOW_LEN; i++) {
            int buffer_i = buffer_index_plus_fftLib(signal_buffer_next_i, i, WINDOW_LEN);
            windowed_signal_in[i] = (double)signal_buffer[buffer_i];
        }

        // Applica la finestra di Hann
        apply_hann_window(windowed_signal_in, windowed_signal_out, WINDOW_LEN);

        // Copia nel vettore complesso per la FFT
        for (int i = 0; i < WINDOW_LEN; i++) {
            fft_input[i].real = windowed_signal_out[i];
            fft_input[i].imag = 0.0;
        }

        
        zero_pad(fft_input, fft_input_padded, WINDOW_LEN, N_PAD);

        // Perform the FFT
        FFT_fftLib(fft_input_padded, N_PAD, 1.0); //NB: If I perform the padding, the index are different

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
            //Without doing the sqrt()
            fft_magnitude[i] = fft_input_padded[i].real * fft_input_padded[i].real + fft_input_padded[i].imag * fft_input_padded[i].imag;

            //Write on the file
#ifdef DUMP_FILE
            if (fft2File)
            {
                fprintf(fft2File, "%f, %f\n", (float)i , fft_magnitude[i]);
                //fprintf(fft2File, "%f, %f, %f\n", (float)i , fft_input_padded[i].real, fft_input_padded[i].imag); 
            }
#endif
        }
        //Finding the max peak
        for (int i = MIN_FREQ_FFT_I + 1; i < MAX_FREQ_FFT_I - 1; i++) {
            if (fft_magnitude[i] > fft_magnitude[i - 1] && fft_magnitude[i] > fft_magnitude[i + 1]) {
                if (fft_magnitude[i] > max_fft_magnitude) {
                    max_fft_magnitude = fft_magnitude[i];
                    dominant_freq_index = i;
                }
            }
        }

        // Calculate the dominant frequency in Hz
        double dominant_freq = (double)dominant_freq_index * SAMPLING_FREQ / N_PAD;

        // Calculate the number of steps based on the dominant frequency
        HR  = dominant_freq * 60;

#ifdef DUMP_FILE
        if (fft2File)
        {
            fflush(fft2File);
            fclose(fft2File);
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

// --- LINEAR INTERPOLATION FUNCTION ---

int fft2_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz){
    int final_HR;
    int HR_results[3] = {-1};
    for(int i=0;i<3;i++){
        HR_results[i] = -1;
    }
    if(delta_ms < 60){
        HR_results[0] = main_algorithm_fft2(delta_ms,ppg,accx,accy,accz);
        final_HR = HR_results[0];
    }
    else if ((delta_ms > 60) & (delta_ms < 100)){
        int delta_ms_interp = delta_ms / 2;
        int ppg_interp = (ppg + ppg_last) / 2;
        int accx_interp = (accx + accx_last) / 2;
        int accy_interp = (accy + accy_last) / 2;
        int accz_interp = (accz + accz_last) / 2;

        HR_results[0] = main_algorithm_fft2(delta_ms_interp, ppg_interp, accx_interp, accy_interp, accz_interp);
        HR_results[1] = main_algorithm_fft2(delta_ms,ppg,accx,accy,accz);
        final_HR = HR_results[1];
    }
    else{
        int delta_ms_interp1 = delta_ms / 3;
        int ppg_interp1 = ppg/3 + ppg_last*2/3;
        int accx_interp1 = accx/3 + accx_last*2/3;
        int accy_interp1 = accy/3 + accy_last*2/3;
        int accz_interp1 = accz/3 + accz_last*2/3;
        int delta_ms_interp2 = delta_ms*2/3;
        int ppg_interp2 = ppg*2/3 + ppg_last/3;
        int accx_interp2 = accx*2/3 + accx_last/3;
        int accy_interp2 = accy*2/3 + accy_last/3;
        int accz_interp2 = accz*2/3 + accz_last/3;

        HR_results[0] = main_algorithm_fft2(delta_ms_interp1, ppg_interp1, accx_interp1, accy_interp1, accz_interp1);
        HR_results[1] = main_algorithm_fft2(delta_ms_interp2, ppg_interp2, accx_interp2, accy_interp2, accz_interp2);
        HR_results[2] = main_algorithm_fft2(delta_ms,ppg,accx,accy,accz);
        final_HR = HR_results[2];
    }
    ppg_last = ppg;
    accx_last = accx;
    accy_last = accy;
    accz_last = accz;

    return final_HR;


}