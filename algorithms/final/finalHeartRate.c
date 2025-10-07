/* ----------------------------------------------------
* FINAL ALGORITHM
* ----------------------------------------------------
* This folder is meant to contain the final algorithm that will be used on the Bangle.js
* 
* Description: The algorithm is based on an NLMS adaptive filter for Motion Artifacts removal, on spectral analysis with FFT and a Kalman filtering 
* for peak tracking.
* The following steps are followed:
* 1. Preprocessing: Band-pass filtering (0.6 - 3.5 Hz) and Standardization (z-score)
* 2. Adaptive filtering with NLMS algorithm, using the three axes of the accelerometer as reference inputs
* 3. State Machine to adapt the algorithm to the current state (rest, motion, transition)
* 4. Spectral analysis with FFT (with zero padding). Peak selection is done exploiting a prior knowledge from the ACC (acceleroeter signal)
* 5. Computation of a confidence index for the selected peak
* 6. HR estimation with a Kalman filter, using the confidence index to update state and measure variance
*/


// --- LIBRARIES ---
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

#include "finalHeartRate.h"


// --- PARAMETERS DEFINITION ---

#define M_PI 3.14159265358979323846    // pi - constant

// For PPG and ACC signal processing: 
#define WINDOW_LEN 128                 // Length of the sliding window to process the signal online
#define WINDOW_STEP 64                 // Overlap: 50% 
#define SAMPLING_FREQ 25               // Sampling frequency (Hz)
#define N_PAD 1024                     // Length of the window after padding (for FFT)

// Heart rate ranges:
// Physiological range: (0.6 - 3.5 Hz) -> (36 - 210 bpm)
#define LOW_FREQ_PHY_I 25              // corresponding to 36 bpm => (freq * N_PAD)/SAMPLING_FREQ = 25: with N_PAD = 1024, freq = 0.61Hz (or 36.6 bpm)
#define HIGH_FREQ_PHY_I 143            // corresponing to 210 bpm => (freq * N_PAD)/SAMPLING_FREQ = 143: with N_PAD = 1024, freq = 3.49.Hz (or 209.5 bpm)
//HR_range (in indexes)
#define DELTA_HR_RANGE 48 //34              // corresponding to 50 bpm => (freq * N_PAD)/SAMPLING_FREQ = 34 with N_PAD = 1024, freq = 0.83Hz (or 49.8 bpm)

// State Machine
#define THRESHOLD_MOTION 300.0
#define STATE_STATIONARY 0
#define STATE_MOTION 1
#define STATE_TRANSITION 2

// For NLMS Adaptive filtering 
#define MU 0.2                         // Learning rate of the NLMS filter

// Confidence coefficient
#define DELTA_PEAK_i 16                // corresponding to 24 bpm or 0.4 Hz => (freq * N_PAD)/SAMPLING_FREQ = 16 with N_PAD = 1024       
#define DELTA_C3_i 8                   // corresponding to 12 bpm or 0.2 Hz => (freq * N_PAD)/SAMPLING_FREQ = 8 with N_PAD = 1024
#define SCALE_C1 10.0                  // Scale factpr Sigmoid c1
#define SCALE_C2 10.0                  // Scale factpr Sigmoid c2
#define MIDPOINT_C1 0.35               // Midpoint Sigmoid c1
#define MIDPOINT_C2 0.65               // Midpoint Sigmoid c2
#define ALPHA 0.45                     // Weight for c1
#define BETA  0.45                     // Weight for c2
#define GAMMA 0.10                     // Weight for c3

// Saving the last spectral peaks
#define N_LAST_SAVED 3

#define EPSILON 1e-10f

// Evaluation of high confidence level
#define THRESHOLD_HIGH_CONFIDENCE 0.87

// To save .csv files 
// #define DUMP_FILE 
#ifdef DUMP_FILE
static int fft_passes = 0;             // counter of how many times the autocorr has been called
#define DUMP_FFT_ST_FILE_NAME "bangle.js-HR/Dump_files/fft_ST"                       // To save the FFT at each window
static FILE *fftSTFile;
#define DUMP_VALUESWIND_ST_FILE_NAME "bangle.js-HR/Dump_files/values_window_C.csv"   // To save some values for each window
static FILE *values_win_C_File;
#define DUMP_SIGNALS_ST_FILE_NAME "bangle.js-HR/Dump_files/list_signals_C.csv"       // To save the overall signals
static FILE *list_signals_C_File;
#define DUMP_RESULTS_ST_FILE_NAME "bangle.js-HR/Dump_files/list_results_C.csv"       // To save the results (Heart rate)
static FILE *list_results_C_File;
#endif


// --- STATIC VARIABLES ---

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

// Band Pass Filter
static BPFilter bpFilter_ppg;
static BPFilter bpFilter_accx;
static BPFilter bpFilter_accy;
static BPFilter bpFilter_accz;

// Standardization
static Stats stats_ppg;
static Stats stats_accx;
static Stats stats_accy;
static Stats stats_accz;
static Stats stats_acc;
static Stats stats_acc_raw;

// NLMS Adaptive Filter
static AdaptFilter_multi NLMS_filter;

// Taking track of the frst window
static bool first_window;
static bool second_window;

// Buffers and counters
static double signal_PPG_buffer[WINDOW_LEN] = {0.0};                                 // Buffer for the standardized PPG signal
static double signal_NLMS_buffer[WINDOW_LEN] = {0.0};                                // Buffer for the output of NLMS filter
static int signal_buffer_next_i = 0;                                                 // Buffer index
static int samples_since_last_HR = 0;                                                // Counter of samples since last HR computation
static accel_t acc_magnitude[WINDOW_LEN] = {0};                                      // Buffer for the magnitude of the acceleration signal
double windowed_signal_in[WINDOW_LEN];                                               // For the fft and zero-padding
double windowed_signal_out[WINDOW_LEN]; 
static complex_number fft_input[WINDOW_LEN];
static complex_number fft_input_padded[N_PAD];

// State Machine
static int state = 0;
static int counter_over = 0;
static int counter_under = 0;

// Model for HR from acceleration prior knowledge (in bpm)
static double a = 0.006;
static double b = 48.15;

// Array with the last spectral peaks and confidences 
static int last_peaks_i[N_LAST_SAVED] = {0};
static double last_peaks[N_LAST_SAVED] = {0.0};
static double last_confs[N_LAST_SAVED] = {0.0};

// Kalman filter
static double HR_freq_est2 = 1.0;      // Kalman Filter 1         
static double P = 0.2;
static double HR_freq_est3 = 1.0;      // Kalman filter 2
static double P_2m = 0.1;

// To have evaluate the confidence value of each measure
static bool high_confidence = false;

//Final value of HR
static int HR = 0;


// --- SOME FUNCTIONS ---

// Model for HR from acceleration prior knowledge (in bpm)
double HR_fom_ACC_double2(double rms, double a, double b){
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

// Update the HR model from the acceleration
// The update is done only on the b parameter, using a gradient update rule
void update_HR_model_from_ACC(double HR_PPG, double HR_ACC, double* a, double* b){
    double learning_rate = 0.25;
    // Keep a unchanged
    *b = *b + learning_rate * (HR_PPG - HR_ACC);
    if (*b < -20.0) *b = -20.0;
    if (*b > 50.0) *b = 50.0;

}

// Funzione sigmoid
double sigmoid_conf2(double x, double midpoint, double scale) {
    return 1.0 / (1.0 + exp(-scale * (x - midpoint)));
}

// Weighted average of the last peaks with their confidence
void weighted_average_with_confidence2(double* peaks, double* confidences, int len, double* weighted_avg, double* C_est) {
    double sum_weighted_peaks = 0.0;
    double sum_conf = 0.0;
    double sum_conf_squared = 0.0;

    for (int i = 0; i < len; i++) {
        sum_weighted_peaks += peaks[i] * confidences[i];
        sum_conf += confidences[i];
        sum_conf_squared += confidences[i] * confidences[i];
    }

    sum_conf += EPSILON;

    *weighted_avg = sum_weighted_peaks / sum_conf;
    *C_est = sum_conf_squared / sum_conf;
}


// --- INIT FUNCTION --- 
void final_heartrate_init()
{
    // Linear interpolation
    ppg_last=0;
    accx_last=0;
    accy_last=0;
    accz_last=0;    

    // Initialization BandPass filters
    BPFilter_init(&bpFilter_ppg);
    BPFilter_init(&bpFilter_accx);
    BPFilter_init(&bpFilter_accy);
    BPFilter_init(&bpFilter_accz);

    // Initialization Standardization
    rolling_stats_reset(&stats_ppg);
    rolling_stats_reset(&stats_accx);
    rolling_stats_reset(&stats_accy);
    rolling_stats_reset(&stats_accz);
    rolling_stats_reset(&stats_acc);
    rolling_stats_reset(&stats_acc_raw);

    // Initialization NLMS adaptive filter
    AdaptFilter_multi_init(&NLMS_filter);


    // Initialize the Buffer to zeros
    for (int i = 0; i < WINDOW_LEN; i++)
    {
        signal_PPG_buffer[i] = 0.0;
        signal_NLMS_buffer[i] = 0.0;
        windowed_signal_in[i] = 0.0;
        windowed_signal_out[i] = 0.0;
    }
    signal_buffer_next_i = 0;

    // Taking track of the first window
    first_window = true;
    second_window = true;

    samples_since_last_HR = 0;
    
    // State Machine 
    state = STATE_STATIONARY;                          // Initialized to STATE_STATIONARY
    counter_over = 0;
    counter_under = 0;

    // Model for HR from acceleration prior knowledge (in bpm)
    a = 0.006;
    b = 48.15;

    // Initialization of the last peaks and confidences
    for (int i=0; i<N_LAST_SAVED;i++){
        last_peaks_i[i] = 0;
        last_peaks[i] = 0.0;
        last_confs[i] = 0.0;
    }

    // Initialization Kalman filters
    kalman_HR_init(&HR_freq_est2, &P);
    kalman_HR_2meas_init(&HR_freq_est3, &P_2m);

    // To evaluate the confidence level of each mesure
    high_confidence = false;

    // Final results initialized to zero
    HR = 0;
    
    // Initialization of the dumping files - open files and writing the headers
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
        fprintf(list_results_C_File, "%s, %s, %s, %s\n", "HR_kalman", "HR_kalman_SF", "a", "b");
    }
#endif

}


// --- MAIN ALGORITHM FUNCTION ---

// The function takes as input the samples already corrected for the missing values.
// As output it gives the estimation of the heart rate on successives windows.

int main_algorithm(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz){

    // Acceleration magnitude - raw
    accel_t acc_raw = (accel_t)sqrt(accx*accx + accy*accy + accz*accz);

    // Applying the Bandpass filters
    BPFilter_put(&bpFilter_ppg, ppg);
    ppg_t ppg_filtered = BPFilter_get(&bpFilter_ppg);
    BPFilter_put(&bpFilter_accx, accx);
    accel_t accx_filtered = BPFilter_get(&bpFilter_accx);
    BPFilter_put(&bpFilter_accy, accy);
    accel_t accy_filtered = BPFilter_get(&bpFilter_accy);
    BPFilter_put(&bpFilter_accz, accz);
    accel_t accz_filtered = BPFilter_get(&bpFilter_accz);

    // Acceleration magnitude - filtered
    accel_t acc_filtered = (accel_t)sqrt(accx_filtered*accx_filtered + accy_filtered*accy_filtered + accz_filtered*accz_filtered);

    // Standardization z-score

    // PPG signal
    rolling_stats_addValue((double)ppg_filtered, &stats_ppg);
    double mean_ppg = rolling_stats_get_mean(&stats_ppg);
    if (rolling_stats_get_variance(&stats_ppg) < 0.0){
        rolling_stats_reset(&stats_ppg);
    }
    double std_ppg = rolling_stats_get_standard_deviation(&stats_ppg);
    if (std_ppg == 0.0) {
        std_ppg = 1.0; 
    }
    double ppg_standardized = ((double)ppg_filtered - mean_ppg) / std_ppg;

    // ACC magnitude signal
    rolling_stats_addValue((double)acc_filtered, &stats_acc);
    double mean_acc = rolling_stats_get_mean(&stats_acc);
    if (rolling_stats_get_variance(&stats_acc) < 0.0){
        rolling_stats_reset(&stats_acc);
    }
    double std_acc = rolling_stats_get_standard_deviation(&stats_acc);
    if (std_acc == 0.0) {
        std_acc = 1.0; 
    }
    double acc_standardized = (double)(acc_filtered - mean_acc) / std_acc;

    // ACCx signal
    rolling_stats_addValue((double)accx_filtered, &stats_accx);
    double mean_accx = rolling_stats_get_mean(&stats_accx);
    if (rolling_stats_get_variance(&stats_accx) < 0.0){
        rolling_stats_reset(&stats_accx);
    }
    double accx_standardized = (double)(accx_filtered - mean_accx) / std_acc;

    // ACCy signa
    rolling_stats_addValue((double)accy_filtered, &stats_accy);
    double mean_accy = rolling_stats_get_mean(&stats_accy);
    if (rolling_stats_get_variance(&stats_accy) < 0.0){
        rolling_stats_reset(&stats_accy);
    }
    double accy_standardized = (double)(accy_filtered - mean_accy) / std_acc;

    // ACCz signal
    rolling_stats_addValue((double)accz_filtered, &stats_accz);
    double mean_accz = rolling_stats_get_mean(&stats_accz);
    if (rolling_stats_get_variance(&stats_accz) < 0.0){
        rolling_stats_reset(&stats_accz);
    }
    double accz_standardized = (double)(accz_filtered - mean_accz) / std_acc;

    // ACC_raw (like this for now)
    rolling_stats_addValue((double)acc_raw, &stats_acc_raw);
    double std_acc_raw = rolling_stats_get_standard_deviation(&stats_acc_raw);
    if (rolling_stats_get_variance(&stats_acc_raw) < 0.0){
        rolling_stats_reset(&stats_acc_raw);
    }
    if (std_acc_raw == 0.0) {
        std_acc_raw = 1.0; 
    }

    //TO CHECK UP
    //Provo così: ATTENZIONE HO AGGIUNTO QUESTA! Ma non è questo il problema!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    std_acc_raw = std_acc;


    /// NLMS Adaptive filtering for MAs removal
    AdaptFilter_multi_put(&NLMS_filter, accx_standardized, accy_standardized, accz_standardized, ppg_standardized);
    double out_NLMS;
    out_NLMS = AdaptFilter_multi_get(&NLMS_filter, MU);

    // Dumping of the signals file 
#ifdef DUMP_FILE
    if (list_signals_C_File)
    {
        if (!fprintf(list_signals_C_File, "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %.13f, %.13f, %.13f, %.13f, %.13f, %.13f\n", ppg, accx, accy, accx, acc_raw, ppg_filtered, accx_filtered, accy_filtered, accz_filtered, acc_filtered, ppg_standardized, accx_standardized, accy_standardized, accz_standardized, out_NLMS, std_acc_raw ))
            puts("error writing file");
        fflush(list_signals_C_File);
    }
#endif

    // Add the ppg and the output of the NLMS to the circular buffer
    signal_PPG_buffer[signal_buffer_next_i] = ppg_standardized;
    signal_NLMS_buffer[signal_buffer_next_i] = out_NLMS;
    acc_magnitude[signal_buffer_next_i] = acc_raw;
    signal_buffer_next_i = (signal_buffer_next_i + 1) % WINDOW_LEN;

    samples_since_last_HR++;


    // After WINDOW_STEP samples, check if it's time to perform the analysis
    if ((samples_since_last_HR == WINDOW_STEP) && (!first_window)){

        // File dump: open the file to plot the fft
#ifdef DUMP_FILE
        fft_passes++;
        char fftSTFileName[100] = DUMP_FFT_ST_FILE_NAME;
        char idxstr[5];
        sprintf(idxstr, "%d", fft_passes);
        strcat(fftSTFileName, idxstr);
        strcat(fftSTFileName, ".csv");
        fftSTFile = fopen(fftSTFileName, "w+");
#endif


        // State Machine

        if (state == STATE_STATIONARY) {
            if (std_acc_raw >= THRESHOLD_MOTION) {
                state = STATE_TRANSITION;
                counter_over = 1;
            } else {
                counter_over = 0; 
            }

        } else if (state == STATE_MOTION) {
            if (std_acc_raw < THRESHOLD_MOTION) {
                state = STATE_TRANSITION;
                counter_under = 1;
            } else {
                counter_under = 0; 
            }

        } else if (state == STATE_TRANSITION) {
            // Transition from STATIONARY to MOTION
            if (counter_over > 0) {
                if (std_acc_raw >= THRESHOLD_MOTION) {
                    counter_over++;
                    if (counter_over >= 4) {
                        state = STATE_MOTION;
                        counter_over = 0;
                    }
                } else {
                    state = STATE_STATIONARY;
                    counter_over = 0;
                }
            }

            // Transition from MOTION to STATIONARY
            else if (counter_under > 0) {
                if (std_acc_raw < THRESHOLD_MOTION) {
                    counter_under++;
                    if (counter_under >= 2) { 
                        state = STATE_STATIONARY;
                        counter_under = 0;
                    }
                } else {
                    state = STATE_MOTION;
                    counter_under = 0;
                }
            }
        }

        // Estimation of the HR from the ACC
        
        //Compute the rms_acc value in double
        double rms_double;
        double cumulative_double = 0;
        for (int i = 0; i < WINDOW_LEN; i++) {
            double val_double = (double)acc_magnitude[i];
            cumulative_double += (val_double * val_double);
        }
        double mean_double = cumulative_double / WINDOW_LEN;
        rms_double = sqrt(mean_double);

        //HR-range with the acceleration prior knowledge
        double HR_ACC = HR_fom_ACC_double2(rms_double, a, b);
        double HR_range_centre_i_double = (HR_ACC*N_PAD) / (SAMPLING_FREQ*60.0);
        double low_freq_HRrange_i_double = HR_range_centre_i_double - DELTA_HR_RANGE/2;
        double high_freq_HRrange_i_double = HR_range_centre_i_double + DELTA_HR_RANGE/2;
       

        samples_since_last_HR = 0;


        // Writing on the dumping file the values of the window
        #ifdef DUMP_FILE
if (values_win_C_File)
    {
        if (!fprintf(values_win_C_File, "%d, %f, %f, %f, %f\n", state, rms_double, HR_range_centre_i_double, low_freq_HRrange_i_double, high_freq_HRrange_i_double )) 
            puts("error writing file");
        fflush(values_win_C_File);
    }
#endif

        // --- PerformingFFT ---
        // Preparing the buffer for the fft
        for (int i = 0; i < WINDOW_LEN; i++) {
            int buffer_i = buffer_index_plus_fftLib(signal_buffer_next_i, i, WINDOW_LEN);
            if (state == STATE_STATIONARY){
                windowed_signal_in[i] = signal_PPG_buffer[buffer_i];
            } else {
                windowed_signal_in[i] = signal_NLMS_buffer[buffer_i];
            }
        }

        // Apply Hann window
        apply_hann_window_fftLib(windowed_signal_in, windowed_signal_out, WINDOW_LEN);

        // Copy on the complex vector for the FFT
        for (int i = 0; i < WINDOW_LEN; i++) {
            fft_input[i].real = windowed_signal_out[i];
            fft_input[i].imag = 0.0;
        }

        // Applying zero-padding
        zero_pad_fftLib(fft_input, fft_input_padded, WINDOW_LEN, N_PAD);

        // Perform the FFT
        FFT_fftLib(fft_input_padded, N_PAD, 1.0); 


        // --- Peak selection ---
        double max_fft_magnitude = 0.0;
        int dominant_freq_index = 0;

        double fft_magnitude[HIGH_FREQ_PHY_I] = {0.0};                    //The first ones - until MIN_FREQ_FFT_I - are initialized to zero and they keep being zero
        double total_power = 0;
        for (int i = LOW_FREQ_PHY_I; i < (HIGH_FREQ_PHY_I+1); i++)
        {
            fft_magnitude[i] = (fft_input_padded[i].real * fft_input_padded[i].real + fft_input_padded[i].imag * fft_input_padded[i].imag) / N_PAD;
            total_power += fft_magnitude[i];

            //Writing on the file the FFT on the physiological range
#ifdef DUMP_FILE
            if (fftSTFile)
            {
                fprintf(fftSTFile, "%f, %f\n", (float)i , fft_magnitude[i]); 
            }
#endif
        }

        // Defining the lower and the higher limit for the peak search
        int lower_limit = LOW_FREQ_PHY_I;
        int higher_limit = HIGH_FREQ_PHY_I;
        if(state == STATE_MOTION){
            lower_limit = (int)ceil(low_freq_HRrange_i_double); 
            higher_limit = (int)floor(high_freq_HRrange_i_double);
            // Clamp to the physiological range
            if (lower_limit <= LOW_FREQ_PHY_I) lower_limit = LOW_FREQ_PHY_I;
            if (higher_limit >= HIGH_FREQ_PHY_I) higher_limit = HIGH_FREQ_PHY_I;

        }

        //Finding the max peak
        for (int i = lower_limit + 1; i < higher_limit; i++) {            //Poi rimettere higher_limit - 1
            if (fft_magnitude[i] > fft_magnitude[i - 1] && fft_magnitude[i] > fft_magnitude[i + 1]) {
                if (fft_magnitude[i] > max_fft_magnitude) {
                    max_fft_magnitude = fft_magnitude[i];
                    dominant_freq_index = i;
                }
            }
        }

        // Check if one peak is found
        bool peak_found = true;
        if (dominant_freq_index == 0.0){
            dominant_freq_index = last_peaks_i[0];                        // If no peak is found, I take the last one NB: ADD SOMETHING IN CASE OF THE FIRST WINDOW
            peak_found = false;
        }

        //Saving the dominant peak in the vector of the last spectral peaks
        for(int i=N_LAST_SAVED-1; i > 0; i--){
            last_peaks_i[i] = last_peaks_i[i-1];
            last_peaks[i] = last_peaks[i-1];
        }
        last_peaks_i[0] = dominant_freq_index;                                                       //Index of the FFT
        last_peaks[0] = (double)(dominant_freq_index*SAMPLING_FREQ / (double)N_PAD);                 //In Hz



        // --- Confidence coefficient ---

        // C1

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

        // C2
        double peak_power2 = 0;
        double total_power_HRrange = 0;
        double total_power_50bpm = 0;
        if (state == STATE_MOTION){ 
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
        if (state == STATE_MOTION){
            if(total_power_HRrange != 0.0) c2 = (double)(peak_power2 / total_power_HRrange);
            else c2 = 0.0f;
        } else {
            if(total_power_50bpm != 0.0) c2 = (double)(peak_power2 / total_power_50bpm);
            else c2 = 0.0f;
        }

        // C3
        double c3 = 0.0;
        if (abs(last_peaks_i[0] - last_peaks_i[1]) <= DELTA_C3_i/2){
            c3 = 0.75;
            if (abs(last_peaks_i[0] - last_peaks_i[2]) <= DELTA_C3_i/2){
                c3 = 1.0;
            }
        }

        // Through sigmoid functions:
        double c1_sigmoid = sigmoid_conf2(c1, MIDPOINT_C1, SCALE_C1);
        double c2_sigmoid = sigmoid_conf2(c2, MIDPOINT_C2, SCALE_C2);

        // Weighted combination
        double coeff_sigmoid = ALPHA * c1_sigmoid + BETA * c2_sigmoid + GAMMA * c3;

        // Penalization for transiztion state
        if (state == STATE_TRANSITION){
            coeff_sigmoid = coeff_sigmoid - 0.4;
        }
        // Clamp between 0 e 1
        if (coeff_sigmoid < 0.0) coeff_sigmoid = 0.0;
        if (coeff_sigmoid > 1.0) coeff_sigmoid = 1.0;

        // If the peak is not found, the cinfidence is zero         //NB: MAYBE CHANGE THIS, I COULD SIMPLY REDUCE THE PREVIOUS CONFIDENCE MAYBE!! NAHH MAYBE IT MAKES SENSE
        if (!peak_found){
            coeff_sigmoid = 0.0;
        }
        

        // --- Spectral tracking ---
        double HR_meas = last_peaks[0];

        //Saving the confidence coefficients in the vector of the last confidences
        for(int i=N_LAST_SAVED-1; i > 0; i--){
            last_confs[i] = last_confs[i-1];
        }
        last_confs[0] = coeff_sigmoid;

        // Kalman Filter 1 
        kalman_HR_estimation(HR_meas, rms_double, coeff_sigmoid, &HR_freq_est2, &P);

        // Kalman Filter 2 (Kalman with two measures: SENSOR FUSION) 
        kalman_HR_estimation_2measures(HR_meas, (HR_ACC/60.0), rms_double, coeff_sigmoid, &HR_freq_est3, &P_2m);

        //Computation of the dinamic range of confidance:
        // Kalman Filter 1
        double HR_low_est_limit = HR_freq_est2 - 1.96 * sqrt(P);
        double HR_high_est_limit = HR_freq_est2 + 1.96 * sqrt(P);

        // Klman Filter 2
        double HR_low_est_limit_2m = HR_freq_est3 - 1.96 * sqrt(P_2m);
        double HR_high_est_limit_2m = HR_freq_est3 + 1.96 * sqrt(P_2m);

        // Final HR value - Without sensor fusion 
        HR  = (int)(HR_freq_est2 * 600); //in bpm x 10

        double HR_k3 = HR_freq_est3*600;
        double HR_k2 = HR_freq_est2*600;

        // Evaluate the confidence level
        int count_confident_meas = 0;
        for (int i=0; i<N_LAST_SAVED;i++){
            if (last_confs[i] >= THRESHOLD_HIGH_CONFIDENCE){
                count_confident_meas++;
            }
        }
        if (count_confident_meas >= 2){
            high_confidence = true;
        } else {
            high_confidence = false;
        }

        // If the confidence is high => update the HR model from the ACC
        if (high_confidence){
            // Compute the weighted average of the last peaks with their confidence
            double HR_freq_average;
            double c_average;
            weighted_average_with_confidence2(last_peaks, last_confs, N_LAST_SAVED, &HR_freq_average, &c_average);

            // Update the model
            update_HR_model_from_ACC(HR_freq_average*60, HR_ACC, &a, &b);

        }
        


        // File dumping
#ifdef DUMP_FILE
        // Writing on the results file
        if (list_results_C_File)
            {
                if (!fprintf(list_results_C_File, "%f, %f, %f, %f\n", HR_k2, HR_k3, a, b )) 
                    puts("error writing file");
                fflush(list_results_C_File);
            }

        // Writing on the fft file some values of the window
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

        // Closing the fft file
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

// --- LINEAR INTERPOLATION FUNCTION ---

int final_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz){
    int final_HR;
    int HR_results[3] = {-1};
    for(int i=0;i<3;i++){
        HR_results[i] = -1;
    }
    if(delta_ms < 60){
        HR_results[0] = main_algorithm(delta_ms,ppg,accx,accy,accz);
        final_HR = HR_results[0];
    }
    else if ((delta_ms > 60) & (delta_ms < 100)){
        int delta_ms_interp = delta_ms / 2;
        int ppg_interp = (ppg + ppg_last) / 2;
        int accx_interp = (accx + accx_last) / 2;
        int accy_interp = (accy + accy_last) / 2;
        int accz_interp = (accz + accz_last) / 2;

        HR_results[0] = main_algorithm(delta_ms_interp, ppg_interp, accx_interp, accy_interp, accz_interp);
        HR_results[1] = main_algorithm(delta_ms,ppg,accx,accy,accz);
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

        HR_results[0] = main_algorithm(delta_ms_interp1, ppg_interp1, accx_interp1, accy_interp1, accz_interp1);
        HR_results[1] = main_algorithm(delta_ms_interp2, ppg_interp2, accx_interp2, accy_interp2, accz_interp2);
        HR_results[2] = main_algorithm(delta_ms,ppg,accx,accy,accz);
        final_HR = HR_results[2];
    }
    ppg_last = ppg;
    accx_last = accx;
    accy_last = accy;
    accz_last = accz;

    return final_HR;


}