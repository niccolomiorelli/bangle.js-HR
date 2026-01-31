/* ----------------------------------------------------
* SPECTRAL TRACKING ALGORITHM
* ----------------------------------------------------
* Description: This was the first version of my final algorithm. In final it is the same, but it is written well.
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
#include "../../utils/adaptive_filter/adaptive_filter.h"
#include "../../utils/kalman_filter/kalman_filter.h"

#include "spectralTrackingHeartRate.h"

#define WINDOW_LEN 128  // sliding window length, better if power of 2 (if we want to switch to FFT), 64 samples = 2.56s , 128 samples = 5.12 s, 256 samples = 10.24s
#define WINDOW_STEP 64 // step of the sliding window, 64 samples = 2.56s
#define SAMPLING_FREQ 25 // sampling frequency of the PPG signal
//PHYSIOLOGICAL RANGE
//In case of padding:
#define LOW_FREQ_PHY_I 25 // 25: with N_PAD = 1024 -> 0.61Hz or 36.6 bpm
#define HIGH_FREQ_PHY_I 143 // 143: with N_PAD = 1024 -> 3.49 Hz or 209.5 bpm
//Delta in index for HR_range
#define DELTA_HR_RANGE 34 //In index  of FFT with zero-padding equal to 1204 -> Corresponding to 50 bpm 

#define M_PI 3.14159265358979323846 // pi
#define N_PAD 1024 //Length of the window after padding

//STATE MACHINE
#define THRESHOLD_MOVIMENTO 300.0f
#define STATO_FERMO 0
#define STATO_MOVIMENTO 1
#define STATO_TRANSIZIONE 2

//for CONFIDENCE
#define DELTA_PEAK_i 16 //Corresponding to 0.4Hz if nfft = 1024
#define DELTA_C3_i 8 //Corresponding to 0.2Hz
//Parameters sigmoid
#define SCALE_C1 10.0
#define SCALE_C2 10.0
#define MIDPOINT_C1 0.35
#define MIDPOINT_C2 0.65
//Weights
#define ALPHA 0.45
#define BETA  0.45
#define GAMMA 0.10

//for saving the last spectral peaks
#define N_LAST_SAVED 5

#define EPSILON 1e-10f

//DUMP FILE: to save .csv files of FFT signals for each window
// #define DUMP_FILE 
#ifdef DUMP_FILE
static int fft_passes = 0; // counter of how many times the autocorr has been called
#define DUMP_FFT_ST_FILE_NAME "bangle.js-HR/Dump_files/fft_ST"
static FILE *fftSTFile;
#define DUMP_VALUESWIND_ST_FILE_NAME "bangle.js-HR/Dump_files/values_window_C.csv"
static FILE *values_win_C_File;
#define DUMP_SIGNALS_ST_FILE_NAME "bangle.js-HR/Dump_files/list_signals_C.csv"
static FILE *list_signals_C_File;
#define DUMP_RESULTS_ST_FILE_NAME "bangle.js-HR/Dump_files/list_results_C.csv"
static FILE *list_results_C_File;
#endif

//For the linear interpolation
typedef struct {
    float ppg;
    float accx;
    float accy;
    float accz;
    uint32_t timestamp_ms; // Tempo assoluto in ms
} Sample;



//for State Machine
static int stato = 0;
static int counter_sopra = 0;
static int counter_sotto = 0;

//Buffers and counters
static double signal_PPG_buffer[WINDOW_LEN] = {0.0};
static double signal_NLMS_buffer[WINDOW_LEN] = {0.0};
static int signal_buffer_next_i = 0;
static int HR = 0; 
static int samples_since_last_HR = 0;

//Acc buffer to compute the rms value
static accel_t acc_magnitude[WINDOW_LEN] = {0};

//Band Pass Filter
static BPFilter bpFilter_ppg;
static BPFilter bpFilter_accx;
static BPFilter bpFilter_accy;
static BPFilter bpFilter_accz;
//Standardization
static Stats stats_ppg;
static Stats stats_accx;
static Stats stats_accy;
static Stats stats_accz;
static Stats stats_acc;
static Stats stats_acc_raw;
//NLMS
static AdaptFilter_multi NLMS_filter;
//Complex number structure
// typedef struct
// {
//     double real;
//     double imag;
// } complex_double;

double windowed_signal_in[WINDOW_LEN]; //For the windowing
double windowed_signal_out[WINDOW_LEN]; 
static complex_number fft_input[WINDOW_LEN];
static complex_number fft_input_padded[N_PAD];

//Array with the last three spectral peaks
static int last_peaks_i[N_LAST_SAVED] = {0};
static double last_peaks[N_LAST_SAVED] = {0.0};
static double last_confs[N_LAST_SAVED] = {0.0};

//Kalman filter
static double HR_freq_est2 = 1.0;
static double P = 0.2;
static double HR_freq_est3 = 1.0;
static double P_2m = 0.1;

//Missing samples
static ppg_t ppg_last;
static accel_t accx_last;
static accel_t accy_last;
static accel_t accz_last;

static bool first_window;
static bool second_window;
// Maybe I have to add it to the fft library, let's see

//Crea la finestra di Hann
void apply_hann_window2(double *input, double *windowed_output, int len) {
    for (int n = 0; n < len; n++) {
        double hann = 0.5 * (1.0 - cos(2.0 * M_PI * n / (len - 1)));
        windowed_output[n] = input[n] * hann;
    }
}


//Padding function: it adds zeros to the input, creating another complex_double vector, so I have the two input ready to be tested
void zero_pad2(complex_number *in, complex_number *out, int N, int N_pad) {
    for (int i = 0; i < N_pad; i++) {
        if (i < N) {
            out[i] = in[i];
        } else {
            out[i].real = 0.0;
            out[i].imag = 0.0;
        }
    }
}

float HR_fom_ACC(accel_t rms){
    float a = 0.0093f;
    float b = 14.96f;
    float out = a*(float)rms + b;
    if (out > 210.0){
        out = 210.0;
    }
    if (out < 36.0)
    {
        out = 36.0;
    }
    return out;
}

double HR_fom_ACC_double(double rms){
    double a = 0.0093;
    double b = 14.96;
    double out = a*rms + b;
    if (out > 210.0){
        out = 210.0;
    }
    if (out < 36.0)
    {
        out = 36.0;
    }
    return out;
}

// Funzione sigmoid
double sigmoid_conf(double x, double midpoint, double scale) {
    return 1.0 / (1.0 + exp(-scale * (x - midpoint)));
}



/// Initialise step counting
void spectralTracking_heartrate_init()
{
    HR = 0;
    samples_since_last_HR = 0;
    // Initialize the signal buffer to zeros
    for (int i = 0; i < WINDOW_LEN; i++)
    {
        signal_PPG_buffer[i] = 0.0;
        signal_NLMS_buffer[i] = 0.0;
        windowed_signal_in[i] = 0.0;
        windowed_signal_out[i] = 0.0;
    }
    signal_buffer_next_i = 0;


    BPFilter_init(&bpFilter_ppg);
    BPFilter_init(&bpFilter_accx);
    BPFilter_init(&bpFilter_accy);
    BPFilter_init(&bpFilter_accz);

    rolling_stats_reset(&stats_ppg);
    rolling_stats_reset(&stats_accx);
    rolling_stats_reset(&stats_accy);
    rolling_stats_reset(&stats_accz);
    rolling_stats_reset(&stats_acc);
    rolling_stats_reset(&stats_acc_raw);

    AdaptFilter_multi_init(&NLMS_filter);

    //State Machine 
    stato = 0; //Initialized to FERMO
    counter_sopra = 0;
    counter_sotto = 0;

    for (int i=0; i<N_LAST_SAVED;i++){
        last_peaks_i[i] = 0;
        last_peaks[i] = 0.0;
        last_confs[i] = 0.0;
    }

    kalman_HR_init(&HR_freq_est2, &P);
    kalman_HR_2meas_init(&HR_freq_est3, &P_2m);

    for(int i =0; i<4;i++){
        ppg_last=0;
        accx_last=0;
        accy_last=0;
        accz_last=0;
    }

    first_window = true;
    second_window = true;
    
#ifdef DUMP_FILE
    fft_passes = 0;
    list_signals_C_File = fopen(DUMP_SIGNALS_ST_FILE_NAME, "w+");
    values_win_C_File = fopen(DUMP_VALUESWIND_ST_FILE_NAME, "w+");
    list_results_C_File = fopen(DUMP_RESULTS_ST_FILE_NAME, "w+");
    //Writing the header of the files
    if(list_signals_C_File){
        fprintf(list_signals_C_File, "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n", "ppg", "accx", "accy", "accx", "acc_raw", "ppg_filtered", "accx_filtered", "accy_filtered", "accz_filtered", "acc_filtered", "ppg_standardized", "accx_standardized", "accy_standardized", "accz_standardized", "out_NLMS", "std_acc_raw");
    }
    if(values_win_C_File){
        fprintf(values_win_C_File, "%s, %s, %s, %s, %s\n", "state" , "acc_rms", "HR_from_ACC", "HR_range_low", "HR_range_high");
    }
    if(list_results_C_File){
        fprintf(list_results_C_File, "%s, %s\n", "HR_kalman", "HR_kalman_SF");
    }
#endif

}

void weighted_average_with_confidence(float* peaks, float* confidences, int len, float* weighted_avg, float* C_est) {
    float sum_weighted_peaks = 0.0f;
    float sum_conf = 0.0f;
    float sum_conf_squared = 0.0f;

    for (int i = 0; i < len; i++) {
        sum_weighted_peaks += peaks[i] * confidences[i];
        sum_conf += confidences[i];
        sum_conf_squared += confidences[i] * confidences[i];
    }

    sum_conf += EPSILON;

    *weighted_avg = sum_weighted_peaks / sum_conf;
    *C_est = sum_conf_squared / sum_conf;
}



int algorithm(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz){


    accel_t acc_raw = (accel_t)sqrt(accx*accx + accy*accy + accz*accz);

 
    //Applying the filter
    BPFilter_put(&bpFilter_ppg, ppg);
    ppg_t ppg_filtered = BPFilter_get(&bpFilter_ppg);
    BPFilter_put(&bpFilter_accx, accx);
    accel_t accx_filtered = BPFilter_get(&bpFilter_accx);
    BPFilter_put(&bpFilter_accy, accy);
    accel_t accy_filtered = BPFilter_get(&bpFilter_accy);
    BPFilter_put(&bpFilter_accz, accz);
    accel_t accz_filtered = BPFilter_get(&bpFilter_accz);

    accel_t acc_filtered = (accel_t)sqrt(accx_filtered*accx_filtered + accy_filtered*accy_filtered + accz_filtered*accz_filtered);

    //Standardization
    //PPG
    rolling_stats_addValue((double)ppg_filtered, &stats_ppg);
    double mean_ppg = rolling_stats_get_mean(&stats_ppg);
    if (rolling_stats_get_variance(&stats_ppg) < 0.0){
        rolling_stats_reset(&stats_ppg);
    }
    double std_ppg = rolling_stats_get_standard_deviation(&stats_ppg);
    if (std_ppg == 0.0) {
        std_ppg = 1.0; // Avoid division by zero
    }
    double ppg_standardized = ((double)ppg_filtered - mean_ppg) / std_ppg;

    //ACC
    rolling_stats_addValue((double)acc_filtered, &stats_acc);
    double mean_acc = rolling_stats_get_mean(&stats_acc);
    if (rolling_stats_get_variance(&stats_acc) < 0.0){
        rolling_stats_reset(&stats_acc);
    }
    double std_acc = rolling_stats_get_standard_deviation(&stats_acc);
    if (std_acc == 0.0) {
        std_acc = 1.0; // Avoid division by zero
    }
    double acc_standardized = (double)(acc_filtered - mean_acc) / std_acc;

    //ACCx
    rolling_stats_addValue((double)accx_filtered, &stats_accx);
    double mean_accx = rolling_stats_get_mean(&stats_accx);
    if (rolling_stats_get_variance(&stats_accx) < 0.0){
        rolling_stats_reset(&stats_accx);
    }
    double accx_standardized = (double)(accx_filtered - mean_accx) / std_acc;

    //ACCy
    rolling_stats_addValue((double)accy_filtered, &stats_accy);
    double mean_accy = rolling_stats_get_mean(&stats_accy);
    if (rolling_stats_get_variance(&stats_accy) < 0.0){
        rolling_stats_reset(&stats_accy);
    }
    double accy_standardized = (double)(accy_filtered - mean_accy) / std_acc;

    //ACCz
    rolling_stats_addValue((double)accz_filtered, &stats_accz);
    double mean_accz = rolling_stats_get_mean(&stats_accz);
    if (rolling_stats_get_variance(&stats_accz) < 0.0){
        rolling_stats_reset(&stats_accz);
    }
    double accz_standardized = (double)(accz_filtered - mean_accz) / std_acc;

    //ACC_raw (like this for now)
    rolling_stats_addValue((double)acc_raw, &stats_acc_raw);
    double std_acc_raw = rolling_stats_get_standard_deviation(&stats_acc_raw);
    if (rolling_stats_get_variance(&stats_acc_raw) < 0.0){
        rolling_stats_reset(&stats_acc_raw);
    }
    if (std_acc_raw == 0.0) {
        std_acc_raw = 1.0; // Avoid division by zero
    }

    //Provo così: ATTENZIONE HO AGGIUNTO QUESTA! Ma non è questo il problema!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    std_acc_raw = std_acc;


    ///////////////////////////////////////////////////////////////////////
    ///NLMS filter
    double mu = 0.2;
    AdaptFilter_multi_put(&NLMS_filter, (double)accx_standardized, (double)accy_standardized, (double)accz_standardized, (double)ppg_standardized);
    double out_NLMS;
    out_NLMS = AdaptFilter_multi_get(&NLMS_filter, mu);

#ifdef DUMP_FILE
    if (list_signals_C_File)
    {
        if (!fprintf(list_signals_C_File, "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %.13f, %.13f, %.13f, %.13f, %.13f, %.13f\n", ppg, accx, accy, accx, acc_raw, ppg_filtered, accx_filtered, accy_filtered, accz_filtered, acc_filtered, ppg_standardized, accx_standardized, accy_standardized, accz_standardized, out_NLMS, std_acc_raw ))
            puts("error writing file");
        fflush(list_signals_C_File);
    }
#endif

    //Add the ppg and the output of the NLMS to the circular buffer
    signal_PPG_buffer[signal_buffer_next_i] = (double)ppg_standardized;
    signal_NLMS_buffer[signal_buffer_next_i] = out_NLMS;
    acc_magnitude[signal_buffer_next_i] = acc_raw;
    signal_buffer_next_i = (signal_buffer_next_i + 1) % WINDOW_LEN;

    samples_since_last_HR++;



    // After WINDOW_STEP samples, check if it's time to perform the analysis
    if ((samples_since_last_HR == WINDOW_STEP) && (!first_window)){

#ifdef DUMP_FILE
        fft_passes++;
        char fftSTFileName[100] = DUMP_FFT_ST_FILE_NAME;
        char idxstr[5];
        sprintf(idxstr, "%d", fft_passes);
        strcat(fftSTFileName, idxstr);
        strcat(fftSTFileName, ".csv");
        fftSTFile = fopen(fftSTFileName, "w+");
#endif


        //STATE_MACHINE:
        if (stato == STATO_FERMO) {
            if (std_acc_raw >= THRESHOLD_MOVIMENTO) {
                stato = STATO_TRANSIZIONE;
                counter_sopra = 1;
            } else {
                counter_sopra = 0; // rimango fermo
            }

        } else if (stato == STATO_MOVIMENTO) {
            if (std_acc_raw < THRESHOLD_MOVIMENTO) {
                stato = STATO_TRANSIZIONE;
                counter_sotto = 1;
            } else {
                counter_sotto = 0; // rimango in movimento
            }

        } else if (stato == STATO_TRANSIZIONE) {
            // Transizione da FERMO a MOVIMENTO
            if (counter_sopra > 0) {
                if (std_acc_raw >= THRESHOLD_MOVIMENTO) {
                    counter_sopra++;
                    if (counter_sopra >= 4) { // 1 + 3 conferme
                        stato = STATO_MOVIMENTO;
                        counter_sopra = 0;
                    }
                } else {
                    stato = STATO_FERMO;
                    counter_sopra = 0;
                }
            }

            // Transizione da MOVIMENTO a FERMO
            else if (counter_sotto > 0) {
                if (std_acc_raw < THRESHOLD_MOVIMENTO) {
                    counter_sotto++;
                    if (counter_sotto >= 2) { // 1 + 1 conferma
                        stato = STATO_FERMO;
                        counter_sotto = 0;
                    }
                } else {
                    stato = STATO_MOVIMENTO;
                    counter_sotto = 0;
                }
            }
        }

        //Regression for estimating HR_range
        //First: compute the rms value
        accel_t rms;
        uint64_t cumulative = 0;
        for (int i = 0; i < WINDOW_LEN; i++) {
            int32_t val = (int32_t) acc_magnitude[i];
            cumulative += (uint64_t)(val * val);
        }
        uint64_t mean = cumulative / WINDOW_LEN;
        rms = (accel_t)sqrt(mean);

        //Compute the rms value in double
        double rms_double;
        double cumulative_double = 0;
        for (int i = 0; i < WINDOW_LEN; i++) {
            double val_double = (double) acc_magnitude[i];
            cumulative_double += (val_double * val_double);
        }
        double mean_double = cumulative_double / WINDOW_LEN;
        rms_double = sqrt(mean_double);



        // int HR_range_centre_i = (int)((HR_fom_ACC(rms)*N_PAD) / (SAMPLING_FREQ*60.0));
        int HR_range_centre_i = (int)round(((HR_fom_ACC(rms)*N_PAD) / (SAMPLING_FREQ*60.0)));
        int low_freq_HRrange_i = HR_range_centre_i - DELTA_HR_RANGE/2;
        int high_freq_HRrange_i = HR_range_centre_i + DELTA_HR_RANGE/2;

        //Doing everything with double
        double HR_range_centre_i_double = (HR_fom_ACC_double(rms_double)*N_PAD) / (SAMPLING_FREQ*60.0);
        double low_freq_HRrange_i_double = HR_range_centre_i_double - DELTA_HR_RANGE/2;
        double high_freq_HRrange_i_double = HR_range_centre_i_double + DELTA_HR_RANGE/2;
        //Now I will try starting from the HR in bpm:
        // double low_freq_HRrange_i_double = HR_range_centre_i_double - (25*1024)/(60*25);
        // double high_freq_HRrange_i_double = HR_range_centre_i_double + (25*1024)/(60*25);
        

        

        samples_since_last_HR = 0;


        #ifdef DUMP_FILE
if (values_win_C_File)
    {
        if (!fprintf(values_win_C_File, "%d, %d, %f, %f, %f\n", stato, rms, HR_range_centre_i_double, low_freq_HRrange_i_double, high_freq_HRrange_i_double )) //Posso plottare sia double che int utile per confronto
            puts("error writing file");
        fflush(values_win_C_File);
    }
#endif

        //Prepare the data for FFT
        for (int i = 0; i < WINDOW_LEN; i++) {
            int buffer_i = buffer_index_plus_fftLib(signal_buffer_next_i, i, WINDOW_LEN);
            if (stato == STATO_FERMO){
                windowed_signal_in[i] = signal_PPG_buffer[buffer_i];
            } else {
                windowed_signal_in[i] = signal_NLMS_buffer[buffer_i];
            }
        }

        // Applica la finestra di Hann
        apply_hann_window2(windowed_signal_in, windowed_signal_out, WINDOW_LEN);

        //Copia nel vettore complesso per la FFT
        for (int i = 0; i < WINDOW_LEN; i++) {
            fft_input[i].real = windowed_signal_out[i];
            fft_input[i].imag = 0.0;
        }

        zero_pad2(fft_input, fft_input_padded, WINDOW_LEN, N_PAD);

        // Perform the FFT
        FFT_fftLib(fft_input_padded, N_PAD, 1.0); //NB: If I perform the padding, the index are different

        // Find the dominant frequency
        double max_fft_magnitude = 0.0;
        int dominant_freq_index = 0;

        //Defining th range 

        double fft_magnitude[HIGH_FREQ_PHY_I] = {0.0}; //The first ones - until MIN_FREQ_FFT_I - are initialized to zero and they keep being zero
        double total_power = 0;
        for (int i = LOW_FREQ_PHY_I; i < (HIGH_FREQ_PHY_I+1); i++)
        {
            fft_magnitude[i] = (fft_input_padded[i].real * fft_input_padded[i].real + fft_input_padded[i].imag * fft_input_padded[i].imag) / N_PAD;
            total_power += fft_magnitude[i];

#ifdef DUMP_FILE
            if (fftSTFile)
            {
                fprintf(fftSTFile, "%f, %f\n", (float)i , fft_magnitude[i]); //Plotto la fft negli indici del range fisiologico
            }
#endif
        }
        //Defining the lower and the higher limit for the peak search
        int lower_limit = LOW_FREQ_PHY_I;
        int higher_limit = HIGH_FREQ_PHY_I;
        if(stato==STATO_MOVIMENTO){
            lower_limit = (int)ceil(low_freq_HRrange_i_double); //NB: This part is important 
            higher_limit = (int)floor(high_freq_HRrange_i_double);
            if (lower_limit <= LOW_FREQ_PHY_I) lower_limit = LOW_FREQ_PHY_I;
            if (higher_limit >= HIGH_FREQ_PHY_I) higher_limit = HIGH_FREQ_PHY_I;

        }
        //Finding the max peak
        for (int i = lower_limit + 1; i < higher_limit; i++) {   //Poi rimettere higher_limit - 1
            if (fft_magnitude[i] > fft_magnitude[i - 1] && fft_magnitude[i] > fft_magnitude[i + 1]) {
                if (fft_magnitude[i] > max_fft_magnitude) {
                    max_fft_magnitude = fft_magnitude[i];
                    dominant_freq_index = i;
                }
            }
        }
        //Checking if I peak is found
        bool peak_found = true;
        if (dominant_freq_index == 0.0){
            dominant_freq_index = last_peaks_i[0];
            peak_found = false;
        }

        //Saving the dominant peak in the vector of the last spectral peaks
        for(int i=N_LAST_SAVED-1; i > 0; i--){
            last_peaks_i[i] = last_peaks_i[i-1];
            last_peaks[i] = last_peaks[i-1];
        }
        last_peaks_i[0] = dominant_freq_index;
        last_peaks[0] = (double)(dominant_freq_index*SAMPLING_FREQ / (double)N_PAD);


        //Computing the CONFIDENCE

        //1. C1
        //I have to compute the power around the peak
        // 1. C1
        double peak_power = 0;
        int peak_lower = dominant_freq_index - DELTA_PEAK_i/2;
        int peak_upper = dominant_freq_index + DELTA_PEAK_i/2;
        //Clamp to valid range
        if (peak_lower < LOW_FREQ_PHY_I) peak_lower = LOW_FREQ_PHY_I;
        if (peak_upper > HIGH_FREQ_PHY_I) peak_upper = HIGH_FREQ_PHY_I;

        for (int i = peak_lower; i <= peak_upper; i++) {
            peak_power += fft_magnitude[i];
        }
        double c1 = (double)(peak_power/total_power);

        // 2. C2
        double peak_power2 = 0;
        double total_power_HRrange = 0;
        double total_power_50bpm = 0;
        if (stato == STATO_MOVIMENTO){ //Non proprio uguale a ciò che viene fatto in python
            int peak_lower2 = dominant_freq_index - DELTA_PEAK_i/2;
            int peak_upper2 = dominant_freq_index + DELTA_PEAK_i/2;
            //Clamp to HR_range
            if (peak_lower2 < lower_limit) peak_lower2 = lower_limit;
            if (peak_upper2 > higher_limit) peak_upper2 = higher_limit;
            for (int i = peak_lower2; i <= peak_upper2; i++) {
                peak_power2 += fft_magnitude[i];
            }

            for (int i = lower_limit; i <= higher_limit; i++) {
                total_power_HRrange += fft_magnitude[i];
            }
        } else {

            int peak_lower2 = dominant_freq_index - DELTA_PEAK_i/2;
            int peak_upper2 = dominant_freq_index + DELTA_PEAK_i/2;
            //Clamp to HR_range
            if (peak_lower2 < LOW_FREQ_PHY_I) peak_lower2 = LOW_FREQ_PHY_I;
            if (peak_upper2 > HIGH_FREQ_PHY_I) peak_upper2 = HIGH_FREQ_PHY_I;
            for (int i = peak_lower2; i <= peak_upper2; i++) {
                peak_power2 += fft_magnitude[i];
            }


            int hr_lower = dominant_freq_index - DELTA_HR_RANGE/2;
            int hr_upper = dominant_freq_index + DELTA_HR_RANGE/2;
            //Clamp to valid range
            if (hr_lower < LOW_FREQ_PHY_I) hr_lower = LOW_FREQ_PHY_I;
            if (hr_upper > HIGH_FREQ_PHY_I) hr_upper = HIGH_FREQ_PHY_I;

            for (int i = hr_lower; i <= hr_upper; i++) {
                total_power_50bpm += fft_magnitude[i];
            }
        }
        double c2;
        if (stato == STATO_MOVIMENTO){
            if(total_power_HRrange != 0.0) c2 = (double)(peak_power2 / total_power_HRrange);
            else c2 = 0.0f;
        } else {
            if(total_power_50bpm != 0.0) c2 = (double)(peak_power2 / total_power_50bpm);
            else c2 = 0.0f;
        }

        //3. C3
        double c3 = 0.0;
        if (abs(last_peaks_i[0] - last_peaks_i[1]) <= DELTA_C3_i/2){
            c3 = 0.75;
            if (abs(last_peaks_i[0] - last_peaks_i[2]) <= DELTA_C3_i/2){
                c3 = 1.0;
            }
        }

        //Through sigmoid functions:
        double c1_sigmoid = sigmoid_conf(c1, MIDPOINT_C1, SCALE_C1);
        double c2_sigmoid = sigmoid_conf(c2, MIDPOINT_C2, SCALE_C2);

        //Weighted combination
        double coeff_sigmoid = ALPHA * c1_sigmoid + BETA * c2_sigmoid + GAMMA * c3;

        if (stato == STATO_TRANSIZIONE){
            coeff_sigmoid = coeff_sigmoid - 0.4;
        }
        //Clip tra 0 e 1
        if (coeff_sigmoid < 0.0) coeff_sigmoid = 0.0;
        if (coeff_sigmoid > 1.0) coeff_sigmoid = 1.0;

        //
        if (!peak_found){
            coeff_sigmoid = 0.0;
        }
        

        //Saving the last confidences (for the moving average)
        for (int i = N_LAST_SAVED - 1; i > 0; i--) {
            last_confs[i]  = last_confs[i - 1];
        }
        last_confs[0] = coeff_sigmoid;



        //FILTERS:
        //a. WEIGHTED MOVING AVERAGE
        float HR_freq_est1;
        float c_est1;
        weighted_average_with_confidence(last_peaks, last_confs, N_LAST_SAVED, &HR_freq_est1, &c_est1);

        //b. Kalman filter
        double HR_meas = HR_freq_est1;
        //Using actual eak measure (I think eventually I will use this)
        HR_meas = last_peaks[0]; 
        double last_confs_mean = 0.0;
        for (int i=0;i<N_LAST_SAVED;i++){
            last_confs_mean += last_confs[i];
        }
        last_confs_mean /= 5.0;
        //oppure se voglio
        last_confs_mean = coeff_sigmoid;
        kalman_HR_estimation(HR_meas, rms_double, last_confs_mean, &HR_freq_est2, &P);

        //c. Kalman filter with two measures: SENSOR FUSION 
        kalman_HR_estimation_2measures(HR_meas, (HR_fom_ACC_double(rms_double)/60.0), rms_double, coeff_sigmoid, &HR_freq_est3, &P_2m);

        //Computation of the dinamic range of confidance:
        //a. Case of Weighted Moving Average
        // float delta_min = 0.05f;  // -> NB: I should put the spectral resolution -> it is 0.024 with zero-padding = 1024. I take twice the sectral resolution just to be sure      
        // float delta_max = 0.83f;  //This is the range of HR_range (I don't know if it makes sense) -> It is 50 bpm  
        // float delta_Hz = delta_min + (1 - c_est1) * (delta_max - delta_min);

        //b. In case of Kalman Filtering
        double HR_low_est_limit = HR_freq_est2 - 1.96 * sqrt(P);
        double HR_high_est_limit = HR_freq_est2 + 1.96 * sqrt(P);

        //c. In case of Kalman filtering with two measures (SENSOR FUSION)
        double HR_low_est_limit_2m = HR_freq_est3 - 1.96 * sqrt(P_2m);
        double HR_high_est_limit_2m = HR_freq_est3 + 1.96 * sqrt(P_2m);

        // Calculate the number of steps based on the dominant frequency
        //Using Kalman filter
        HR  = (int)(HR_freq_est3 * 600); //in bpm x 10
        double HR_k3 = HR_freq_est3*600;
        double HR_k2 = HR_freq_est2*600;
        // HR = dominant_freq_index;
        // HR = (int)(HR_freq_est1*600);
        // HR = (int)(last_peaks[0]*600);

        #ifdef DUMP_FILE
if (list_results_C_File)
    {
        if (!fprintf(list_results_C_File, "%f, %f\n", HR_k2, HR_k3 )) 
            puts("error writing file");
        fflush(list_results_C_File);
    }
#endif


        ///////////////////////////////////////////

#ifdef DUMP_FILE
        if (fftSTFile)
        {
            fprintf(fftSTFile, "%s, %d\n", "peak_index", dominant_freq_index);
            fprintf(fftSTFile, "%s, %f\n", "c1" , c1_sigmoid);
            fprintf(fftSTFile, "%s, %f\n", "c2" , c2_sigmoid);
            fprintf(fftSTFile, "%s, %f\n", "c3" , c3);
            fprintf(fftSTFile, "%s, %f\n", "coeff" , coeff_sigmoid);
            fprintf(fftSTFile, "%s, %d\n", "Lower_limit", lower_limit);
            fprintf(fftSTFile, "%s, %d\n", "Higher_limit", higher_limit);
        }
#endif
#ifdef DUMP_FILE
        if (fftSTFile)
        {
            fflush(fftSTFile);
            fclose(fftSTFile);
        }
#endif        




        
    }
    else if ((samples_since_last_HR==WINDOW_STEP) && (first_window)){
        first_window = false;
        samples_since_last_HR=0;
    }

    
    // Return the HR*10
    return HR;

}


int spectralTracking_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz){
    int final_HR;
    int HR_results[3] = {-1};
    for(int i=0;i<3;i++){
        HR_results[i] = -1;
    }
    if(delta_ms < 60){
        HR_results[0] = algorithm(delta_ms,ppg,accx,accy,accz);
        final_HR = HR_results[0];
    }
    else if ((delta_ms > 60) & (delta_ms < 100)){
        int delta_ms_interp = delta_ms / 2;
        int ppg_interp = (ppg + ppg_last) / 2;
        int accx_interp = (accx + accx_last) / 2;
        int accy_interp = (accy + accy_last) / 2;
        int accz_interp = (accz + accz_last) / 2;

        HR_results[0] = algorithm(delta_ms_interp, ppg_interp, accx_interp, accy_interp, accz_interp);
        HR_results[1] = algorithm(delta_ms,ppg,accx,accy,accz);
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

        HR_results[0] = algorithm(delta_ms_interp1, ppg_interp1, accx_interp1, accy_interp1, accz_interp1);
        HR_results[1] = algorithm(delta_ms_interp2, ppg_interp2, accx_interp2, accy_interp2, accz_interp2);
        HR_results[2] = algorithm(delta_ms,ppg,accx,accy,accz);
        final_HR = HR_results[2];
    }
    ppg_last = ppg;
    accx_last = accx;
    accy_last = accy;
    accz_last = accz;

    return final_HR;


}

