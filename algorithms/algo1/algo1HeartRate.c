/* ----------------------------------------------------
* ALGO1 ALGORITHM
* ----------------------------------------------------
* Description:
* The algorithm takes overlapping windows of 5.12s (128 samples @ 25 Hz) with step of 2.56s (64 samples @ 25Hz) and applies the Fast Fourier Transform.
* To change the window length -> WINDOW_LEN. To change the overlapping -> WINDOW_STEP   
* Note: if the WINDOW_len is changed, the MIN_FREQ_FFT_I and MAX_FREQ_FFT_I must be changed too.
* The algorithm searches the peak of the Fourier transform, limiting the search between the indexes [3, 22],
* corresponding to the range of HR [35bpm, 246 bpm] or [0.58Hz, 4.29Hz] -> To change this parameters: MIN_FREQ_FFT_I, MAX_FREQ_FFT_I. Formula: freq = i*f_sample/N
* The algorithm performs a median filter on the last 8 output -> HRM_HIST_LEN and HRM_MEDIAN_LEN
* 
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>


#include "../../utils/bandpass_filter/bandpass_filter.h"
#include "../../utils/adaptive_filter/adaptive_filter.h"
#include "../../types.h"
#include "../../utils/fft_library/fft_library.h"


#define WINDOW_LEN 128  // sliding window length, better if power of 2 (if we want to switch to FFT), 64 samples = 2.56s , 128 samples = 5.12 s, 256 samples = 10.24s
#define WINDOW_STEP 64 // step of the sliding window, 64 samples = 2.56s
#define SAMPLING_FREQ 25 // sampling frequency of the PPG signal
#define MIN_FREQ_FFT_I 4 // index of the FFT corresponding to the minimum heart rate -> 3: corresponds to 35 bpm (if N=256  -> 5 corresponds to 30bpm)
#define MAX_FREQ_FFT_I 15 // index of the FFT corresponding to the maximum heart rate -> 22: correspnds to 246 bpm (if N=256  -> 43 corresponds to 258bpm)
//#define M_PI 3.14159265358979323846 // pi
#define TOP_N 3
// Buffers and counters
static float ppg_buffer[WINDOW_LEN] = {0};
static accel_t acc_buffer[WINDOW_LEN] ={0};
static int buffer_next_i = 0;
static int HR = 0; 
static int samples_since_last_HR = 0;

//DUMP FILE: to save .csv files of FFT signals for each window
//#define DUMP_FILE 
#ifdef DUMP_FILE
static int fft_passes = 0; // counter of how many times the autocorr has been called
#define DUMP_FFT_FILE_NAME "fft"
static FILE *fftFile;
#define DUMP_ALGO1_FILE_NAME "algo1.csv"
static FILE *algo1File;
#endif

//For the median filter
//Buffer per implementare media mobile o mediana sui risultati
#define HRM_HIST_LEN 8
#define HRM_MEDIAN_LEN 4
static int fft_results[HRM_HIST_LEN] = {0};
static int fft_results_index = 0;

//Band Pass Filter
static BPFilter bpFilter_ppg;
static BPFilter bpFilter_acc; //For acc_magnitude
static BPFilter bpFilter_accx;
static BPFilter bpFilter_accy;
static BPFilter bpFilter_accz;

//Adaptive Filter
static AdaptFilter adaptFilter1;
static AdaptFilter adaptFilter2;
static AdaptFilter adaptFilter3;

static complex_number fft_input_ppg[WINDOW_LEN];
static complex_number fft_input_acc[WINDOW_LEN];

// For mean and std
int8_t count  = 0;
double mean_ppg = 0;
double variance_ppg = 0;
double std_ppg = 0;
double mean_acc[4] = {0};
double variance_acc[4] = {0};
double std_acc[4] = {0};
/// Initialise step counting
void algo1_heartrate_init()
{
    HR = 0;
    samples_since_last_HR = 0;
    // Initialize the signal buffer to zeros
    for (int i = 0; i < WINDOW_LEN; i++)
    {
        ppg_buffer[i] = 0.0f;
        acc_buffer[i]=0;
    }
    buffer_next_i = 0;

    //The following buffer contains the HR outputs of the algorithm
    for (int i = 0; i < HRM_HIST_LEN; i++)
    {
        fft_results[i] = 0;;
    }
    fft_results_index = 0;

    BPFilter_init(&bpFilter_ppg);
    BPFilter_init(&bpFilter_acc);
    BPFilter_init(&bpFilter_accx);
    BPFilter_init(&bpFilter_accy);
    BPFilter_init(&bpFilter_accz);

    //AdaptFilter
    AdaptFilter_init(&adaptFilter1);
    AdaptFilter_init(&adaptFilter2);
    AdaptFilter_init(&adaptFilter3);

#ifdef DUMP_FILE
    fft_passes = 0;
    algo1File = fopen(DUMP_ALGO1_FILE_NAME, "w+");
#endif

    count = 0;

}

int algo1_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz){

      
    // Applying the filter on ppg signal
    BPFilter_put(&bpFilter_ppg, ppg);
    ppg_t ppg_filtered = BPFilter_get(&bpFilter_ppg);

    // Applying the filter on acc signal 
    accel_t acc_filtered[4]; //0 -> magnitude, 1-> accx, 2-> accy, 3-> accz
    accel_t acc_magnitude = sqrt(accx*accx + accy*accy + accz*accz);
    BPFilter_put(&bpFilter_acc, acc_magnitude);
    acc_filtered[0] = BPFilter_get(&bpFilter_acc);
    // Applying the filter on accx signal 
    BPFilter_put(&bpFilter_accx, accx);
    acc_filtered[1] = BPFilter_get(&bpFilter_accx);
    // Applying the filter on accy signal 
    BPFilter_put(&bpFilter_accy, accy);
    acc_filtered[2] = BPFilter_get(&bpFilter_accy);
    // Applying the filter on accz signal 
    BPFilter_put(&bpFilter_accz, accz);
    acc_filtered[3] = BPFilter_get(&bpFilter_accz);


    
    // Computing the mean and std
    double alpha = 0.001;
    count++;
    if (count ==1){
        mean_ppg = (double)ppg_filtered;
        variance_ppg = 0.0;
        std_ppg = 0.0;

        for(int i=0;i<4;i++){
            mean_acc[i] = (double)acc_filtered[i];
            variance_acc[i] = 0.0;
            std_acc[i] = 0.0;
        }
        
    }
    else{
        mean_ppg = (alpha*(double)ppg_filtered + (1-alpha)*mean_ppg);
        variance_ppg = (alpha*((double)ppg_filtered - mean_ppg)*((double)ppg_filtered -mean_ppg) + (1-alpha)*variance_ppg);
        std_ppg = (sqrt(variance_ppg));

        for(int i=0;i<4;i++){
            mean_acc[i] = (alpha*(double)acc_filtered[i] + (1-alpha)*mean_acc[i]);
            variance_acc[i] = (alpha*((double)acc_filtered[i] - mean_acc[i])*((double)acc_filtered[i] - mean_acc[i]) + (1-alpha)*variance_acc[i]);
            std_acc[i] = (sqrt(variance_acc[i]));
        }
        
        count=7;
    }
    
    double ppg_standardized;
    double acc_standardized[4];
    //Standardization of the signal:
    if (std_ppg != 0.0){
        ppg_standardized = (((double)ppg_filtered - mean_ppg)/(std_ppg));
    } else {
        ppg_standardized = (((double)ppg_filtered - mean_ppg)/(1.0));
    }
    for(int i =0;i<4;i++){
        if (std_acc[i] != 0.0){
            acc_standardized[i] = (((double)acc_filtered[i] - mean_acc[i])/(std_acc[i]));
        } else {
            acc_standardized[i] = (((double)acc_filtered[i] - mean_acc[i])/(1.0));
        }
    }
    
    /*
    //Standardized offline
    float ppg_standardized;
    float acc_standardized[4];

    ppg_standardized = (((float)ppg_filtered - 1818.3757471380814)/(580.5866095595597));
    acc_standardized[0]=(float)acc_filtered[0];;
    acc_standardized[1] = (((float)acc_filtered[1] + 2470.8776719683924)/(592.5991259010865));
    acc_standardized[2] = (((float)acc_filtered[2] - 1031.5485766386385)/(529.6285970003521));
    acc_standardized[3] = (((float)acc_filtered[3] + 1141.1342315874786)/(567.7782822616927));
    */

    //Adaptive Filter - Three stages
    //1.
    AdaptFilter_put(&adaptFilter1, acc_standardized[1], ppg_standardized);
    double e1 = AdaptFilter_get(&adaptFilter1);
    //2.
    AdaptFilter_put(&adaptFilter2, acc_standardized[2], e1);
    double e2 = AdaptFilter_get(&adaptFilter2);
    //3.
    AdaptFilter_put(&adaptFilter3, acc_standardized[3], e2);
    double ppg_adaptfilt = AdaptFilter_get(&adaptFilter3);



    // Add the 2 signals to the circular buffer
    ppg_buffer[buffer_next_i] = ppg_adaptfilt;
    acc_buffer[buffer_next_i] = acc_filtered[0];
    buffer_next_i = (buffer_next_i + 1) % WINDOW_LEN;

#ifdef DUMP_FILE
        if (algo1File)
        {
            if (!fprintf(algo1File, "%f, %f, %f, %f, %f, %f, %f, %f, %f,%f, %f, %f, %f, %f, %f, %f, %f\n", (double)ppg_filtered, ppg_standardized, mean_ppg, std_ppg, (double)acc_magnitude, (double)acc_filtered[0],mean_acc[0],std_acc[0],ppg_adaptfilt,acc_standardized[0], (double)acc_filtered[1], (double)acc_filtered[2], (double)acc_filtered[3],(double)ppg,(double)accx, (double)accy, (double)accz))
                puts("error writing file");
            fflush(algo1File);
        }
#endif

    samples_since_last_HR++;

    // After WINDOW_STEP samples, check if it's time to perform the analysis
    if (samples_since_last_HR >= WINDOW_STEP)
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

        

        // Prepare the data for FFT (both ppg and acc)
        for (int i = 0; i < WINDOW_LEN; i++)
        {
            int buffer_i = buffer_index_plus_fftLib(buffer_next_i, i, WINDOW_LEN);
            fft_input_ppg[i].real = (double)ppg_buffer[buffer_i];
            fft_input_ppg[i].imag = 0.0;
            fft_input_acc[i].real = (double)acc_buffer[buffer_i];
            fft_input_acc[i].imag = 0.0;
            
        }

        // Perform the FFT
        FFT_fftLib(fft_input_ppg, WINDOW_LEN, 1.0);
        FFT_fftLib(fft_input_acc, WINDOW_LEN, 1.0);

        // Find the dominant frequency
        /*
        double max_fft_magnitude_ppg = 0.0;
        int dominant_freq_index_ppg = 0;
        double max_fft_magnitude_acc = 0.0;
        int dominant_freq_index_acc = 0;
        */


        double fft_magnitude_ppg[MAX_FREQ_FFT_I] = {0.0}; //The first ones - until MIN_FREQ_FFT_I - are initialized to zero and they keep being zero
        double fft_magnitude_acc[MAX_FREQ_FFT_I] = {0.0};
        

        for (int i = MIN_FREQ_FFT_I; i < MAX_FREQ_FFT_I; i++)
        {            
            fft_magnitude_ppg[i] = sqrt(fft_input_ppg[i].real * fft_input_ppg[i].real + fft_input_ppg[i].imag * fft_input_ppg[i].imag);
            fft_magnitude_acc[i] = sqrt(fft_input_acc[i].real * fft_input_acc[i].real + fft_input_acc[i].imag * fft_input_acc[i].imag);

            //Write on the file
#ifdef DUMP_FILE
            if (fftFile)
            {
                fprintf(fftFile, "%d, %f, %f\n", i, fft_magnitude_ppg[i],fft_magnitude_acc[i]);
            }
#endif
        }
        //Finding the two highest peak in the FFT of ACC
        double max_acc[2];
        max_acc[0] = -1.0; //Highest peak
        max_acc[1] = -1.0; //Second highest peak
        int idx_acc[2];
        idx_acc[0] = -1; //Index of the highest peak
        idx_acc[1] = -1; //Index of the second highest peak

        //MIN_FREQ_FFT_I = 4 -> 0.78Hz so 46bpm
        //MAX_FREQ_FFT_I = 14 -> 2.73Hx so 164bpm
        int found_peaks = 0;
        for (int i = MIN_FREQ_FFT_I + 1; i < MAX_FREQ_FFT_I - 1; ++i) {
            if (fft_magnitude_acc[i] > fft_magnitude_acc[i - 1] && fft_magnitude_acc[i] > fft_magnitude_acc[i + 1]) {
                
                //Local peak found
                float current_peak = fft_magnitude_acc[i];

                if (current_peak > max_acc[0]) {
                    max_acc[1] = max_acc[0];
                    idx_acc[1] = idx_acc[0];

                    max_acc[0] = current_peak;
                    idx_acc[0] = i;
                } else if (current_peak > max_acc[1]) {
                    max_acc[1] = current_peak;
                    idx_acc[1] = i;
                }

                found_peaks++;
            }
        }   


        //NB: The resolution with N=128 is around 0.2Hz

        //Searching the peak in the clean PPG signal
        double max_fft_magnitude_ppg = 0.0;
        int dominant_freq_index_ppg = 0;
        int found_peak_ppg = 0;
        for (int i = MIN_FREQ_FFT_I-1; i < MAX_FREQ_FFT_I+1; i++)
        {
            //Making sure that the peak is not in the same position of the two peaks found in the ACC signal (+- 1 index)

            if((i != idx_acc[0]) && (i != idx_acc[1])){
                if ((fft_magnitude_ppg[i] > fft_magnitude_ppg[i-1]) && (fft_magnitude_ppg[i] > fft_magnitude_ppg[i+1]))
                {
                    //Local peak found
                    float current_peak = fft_magnitude_ppg[i];

                    if (current_peak > max_fft_magnitude_ppg)
                    {
                        max_fft_magnitude_ppg = current_peak;
                        dominant_freq_index_ppg = i;
                    }

                    found_peak_ppg++;
                }
            }
        }
        //If no peak has been found, take a fake one @ 70 bpm
        if(found_peak_ppg == 0){
            dominant_freq_index_ppg = 6; //70 bpm, this could happen at the beginning (no peak in the FFT of the ppg signal)
        }

        /*
        //Finding the top 3 peaks in both signals

        int dominant_freq_index_ppg[TOP_N] = {0};
        double top_fft_magnitude_ppg[TOP_N] = {0.0};

        int dominant_freq_index_acc[TOP_N] = {0};
        double top_fft_magnitude_acc[TOP_N] = {0.0};
        

        for (int i = MIN_FREQ_FFT_I; i < MAX_FREQ_FFT_I; i++)
        {
            //PPG
            double mag_ppg = fft_magnitude_ppg[i];
            for (int j = 0; j < TOP_N; j++) {
                if (mag_ppg > top_fft_magnitude_ppg[j]) {
                    for (int k = TOP_N - 1; k > j; k--) {
                        top_fft_magnitude_ppg[k] = top_fft_magnitude_ppg[k - 1];
                        dominant_freq_index_ppg[k] = dominant_freq_index_ppg[k - 1];
                    }
                    top_fft_magnitude_ppg[j] = mag_ppg;
                    dominant_freq_index_ppg[j] = i;
                    break;
                }
            }

            //ACC
            float mag_acc = fft_magnitude_acc[i];
            for (int j = 0; j < TOP_N; j++) {
                if (mag_acc > top_fft_magnitude_acc[j]) {
                    for (int k = TOP_N - 1; k > j; k--) {
                        top_fft_magnitude_acc[k] = top_fft_magnitude_acc[k - 1];
                        dominant_freq_index_acc[k] = dominant_freq_index_acc[k - 1];
                    }
                    top_fft_magnitude_acc[j] = mag_acc;
                    dominant_freq_index_acc[j] = i;
                    break;
                }
            }
        }
        */
        /*
        for (int i = MIN_FREQ_FFT_I; i < MAX_FREQ_FFT_I; i++)
        {
            if (fft_magnitude_ppg[i] > max_fft_magnitude_ppg)
            {
                max_fft_magnitude_ppg = fft_magnitude_ppg[i];
                dominant_freq_index_ppg = i;
            }
            if (fft_magnitude_acc[i] > max_fft_magnitude_acc)
            {
                max_fft_magnitude_acc = fft_magnitude_acc[i];
                dominant_freq_index_acc = i;
            }
        }
        */




        // Calculate the dominant frequency in Hz
        double dominant_freq_ppg = (double)dominant_freq_index_ppg * SAMPLING_FREQ / WINDOW_LEN;
        

        // Calculate the number of steps based on the dominant frequency
        HR  = dominant_freq_ppg * 60;
        

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