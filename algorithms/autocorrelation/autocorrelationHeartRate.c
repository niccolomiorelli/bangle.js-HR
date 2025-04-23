
/* ----------------------------------------------------
* AUTOCORRELATION ALGORITHM
* ----------------------------------------------------
* Description:
* The algorithm takes non overlapping windows of 5s (125 samples @ 25 Hz) and performs the autocorrelation.
* To change the window length -> NUM_TUPLES 
* The algorithm searches the second peak of the autocorrelation, limiting the search between the [8, 50] lags,
* corresponding to the range of HR [30bpm, 180bpm] or [0.5Hz, 3Hz] -> To change this parameters: NUM_AUTOCORR_LAGS, FIRST_AUTOCORR_PEAK_LAG
* The algorithm currently finds the peak with a peak detection on the original autocorrelation signal. 
* The second replica uses a derivative filter to find the peak.
* The algorithm can also perform a logic to verify the detected peak (for now it is commentated/disabled).
* The algorithm inspired by the algorithm in "https://github.com/VirginiaSek/Algo" used for step counting from accelerometer signal.
*
* Modifications have been made in order to adapt the algorithm to the HR computation from PPG signal.
*
*/



#include "autocorrelationHeartRate.h"
#include <stdint.h>
#include <stdio.h>  
#include <string.h>
#include <stdbool.h>

#include "../../utils/bandpass_filter/bandpass_filter.h"
#include "../../types.h"

#define SAMPLING_RATE 25  // 25 hz sampling rate
#define NUM_TUPLES 128 //125    // 5 seconds worth of data
#define WINDOW_LENGTH (NUM_TUPLES / SAMPLING_RATE) // window length in seconds: 5

static int autocorrelation_HR; 
static int autocorrelation_HR_temp = 0;

//Da togliere:
//static ppg_t autocorr_buffer[NUM_TUPLES]; //Buffer where I put the initial signal
//static long autocorr_buffer_index = 0;


#define NUM_AUTOCORR_LAGS  51      // number of lags to calculate for autocorrelation. At a minimum heart rate of 30 bpm -> 0.5 beat/s -> 2 / sampling_period(0.04s) = 2s * 25Hz = 50 lags 
#define FIRST_AUTOCORR_PEAK_LAG 8 // corresponds to the first feasible autocorrelation lag -> at a max heart rate of 3 beats /s -> 0.333s / sampling_period(0.04s) = 8.25 -> 8 lags

// Vrification logic
#define AUTOCORR_DELTA_AMPLITUDE_THRESH 1e4 // this is the min delta between peak and trough of autocorrelation peak
#define AUTOCORR_MIN_HALF_LEN 2             // this is the min number of points the autocorrelation peak should be on either side of the peak

//DUMP FILE: used to print the signal in .csv files
// #define DUMP_FILE 
#ifdef DUMP_FILE
static int autocorr_passes = 0; // counter of how many times the autocorr has been called
#define DUMP_REMOVED_MEAN_FILE_NAME "removed_mean.csv"
#define DUMP_AUTOCORRELATION_FILE_NAME "autocorrelation"
#define DUMP_MAGNITUDE_FILE_NAME "magnitude.csv"
#define DUMP_FILTERED_FILE_NAME "filtered.csv"
static FILE *autocorrelationFile;
static FILE *magnitudeFile;
static FILE *filteredFile;
static FILE *removedMeanFile;
static FILE *autocorrelationFile;
#endif

static ppg_t autocorr_buffer[NUM_TUPLES];
static long autocorr_buffer_index = 0;
static ppg_t filtered_buff[NUM_TUPLES] = {0};    // low pass filtered data
static int64_t autocorr_buff[NUM_AUTOCORR_LAGS] = {0}; // autocorrelation results
// static int64_t deriv[NUM_AUTOCORR_LAGS] = {0};         // derivative

//Buffer to implement median filter on the estimated HR
#define HRM_HIST_LEN 4
#define HRM_MEDIAN_LEN 2
static int autocorr_results[HRM_HIST_LEN] = {0};
static int autocorr_results_index = 0;

//BandPass filter
static BPFilter bpFilter;

void autocorrelation_heartrate_init()
{
    BPFilter_init(&bpFilter);
    //Putting everything to zero
    for (int i = 0; i < NUM_AUTOCORR_LAGS; i++)
    {
        autocorr_buff[i] = 0;
        //deriv[i] = 0;
    }
    autocorrelation_HR = 0;
    for (int i = 0; i < NUM_TUPLES; i++)
    {
        autocorr_buffer[i] = 0;;
    }
    autocorr_buffer_index = 0;
    for (int i = 0; i < HRM_HIST_LEN; i++)
    {
        autocorr_results[i] = 0;;
    }
    autocorr_results_index = 0;


#ifdef DUMP_FILE
    magnitudeFile = fopen(DUMP_MAGNITUDE_FILE_NAME, "w");
    filteredFile = fopen(DUMP_FILTERED_FILE_NAME, "w");
    removedMeanFile = fopen(DUMP_REMOVED_MEAN_FILE_NAME, "w");

    autocorr_passes = 0;
#endif

}   

// Sends the buffer to a band pass filter
static void bandpass_buffer(ppg_t *input_buffer, ppg_t *output_buffer)
{
    for (int i = 0; i < NUM_TUPLES; i++)
    {
        BPFilter_put(&bpFilter, input_buffer[i]);
        output_buffer[i] = BPFilter_get(&bpFilter);

#ifdef DUMP_FILE
        if (filteredFile)
        {
            fprintf(filteredFile, "%d\n", output_buffer[i]);
            fflush(filteredFile);
        }
#endif
    }
}

// Remove mean from filtered data
static void remove_mean(ppg_t *buffer)
{
    int64_t sum_temp = 0; //The sum_temp is added to be sure to convert correctly between different data type (probably it is not necessary)
    ppg_t sum = 0;
    uint16_t i;
    for (i = 0; i < NUM_TUPLES; i++)
    {
        sum_temp += (int64_t)buffer[i];
    }
    sum_temp = (int64_t)((float)(sum_temp) / (float)(NUM_TUPLES));
    sum = (ppg_t)sum_temp;
    for (i = 0; i < NUM_TUPLES; i++)
    {
        buffer[i] -= sum;

#ifdef DUMP_FILE
        if (removedMeanFile)
        {
            fprintf(removedMeanFile, "%d\n", buffer[i]);
        }
#endif
    }

#ifdef DUMP_FILE
    if (removedMeanFile)
        fflush(removedMeanFile);
#endif
}

// Autocorrelation function: it computes the autocorrelatio for lags between 0 and NUM_AUTOCORR_LAGS
static void autocorr(ppg_t *buffer, int64_t *autocorr_buff)
{

#ifdef DUMP_FILE
    autocorr_passes++;
    char autocorrFileName[100] = DUMP_AUTOCORRELATION_FILE_NAME;
    char idxstr[5];
    sprintf(idxstr, "%d", autocorr_passes);
    strcat(autocorrFileName, idxstr);
    strcat(autocorrFileName, ".csv");
    autocorrelationFile = fopen(autocorrFileName, "w+");
#endif


    uint8_t lag;
    uint16_t i;
    int64_t temp_ac;
    for (lag = 0; lag < NUM_AUTOCORR_LAGS; lag++)
    {
        temp_ac = 0;
        for (i = 0; i < NUM_TUPLES - lag; i++)
        {
            temp_ac += (int64_t)buffer[i] * (int64_t)buffer[i + lag];
        }
        autocorr_buff[lag] = temp_ac;

#ifdef DUMP_FILE
        if (autocorrelationFile)
        {
            fprintf(autocorrelationFile, "%u, %lld\n", lag, autocorr_buff[lag]);
        }
#endif
    }
#ifdef DUMP_FILE
    if (autocorrelationFile)
    {
        fflush(autocorrelationFile);
        fclose(autocorrelationFile);
    }
#endif
}

// take a look at the original autocorrelation signal at index i and see if
// it's a real peak or if it's just a fake "noisy" peak corresponding to
// non-walking. Basically just count the number of points of the
// autocorrelation peak to the right and left of the peak. this function gets
// the number of points to the right and left of the peak, as well as the delta amplitude
static void get_autocorr_peak_stats(int64_t *autocorr_buff, uint8_t *neg_slope_count, int64_t *delta_amplitude_right, uint8_t *pos_slope_count, int64_t *delta_amplitude_left, uint8_t peak_ind)
{

    // first look to the right of the peak. walk forward until the slope begins decreasing (increasing in my opinion)
    uint8_t neg_slope_ind = peak_ind;
    uint16_t loop_limit = NUM_AUTOCORR_LAGS - 1;
    while ((autocorr_buff[neg_slope_ind + 1] - autocorr_buff[neg_slope_ind] < 0) && (neg_slope_ind < loop_limit))
    {
        *neg_slope_count = *neg_slope_count + 1;
        neg_slope_ind = neg_slope_ind + 1;
    }

    // get the delta amplitude between peak and right trough (negative peak)
    *delta_amplitude_right = autocorr_buff[peak_ind] - autocorr_buff[neg_slope_ind];

    // next look to the left of the peak. walk backward until the slope begins increasing (decreasing in my opinion)
    uint8_t pos_slope_ind = peak_ind;
    loop_limit = 0;
    while ((autocorr_buff[pos_slope_ind] - autocorr_buff[pos_slope_ind - 1] > 0) && (pos_slope_ind > loop_limit))
    {
        *pos_slope_count = *pos_slope_count + 1;
        pos_slope_ind = pos_slope_ind - 1;
    }

    // get the delta amplitude between the peak and the left trough
    *delta_amplitude_left = autocorr_buff[peak_ind] - autocorr_buff[pos_slope_ind];
}



// Function that returns the period in seconds of the signal (interval) -> it corresponds to the lag of the second peak of autocorrelation divided by f_sample
static float autcorr_count_beats(ppg_t *ppg_buffer)
{

#ifdef DUMP_FILE
    for (int i = 0; i < NUM_TUPLES; i++)
    {
        if (magnitudeFile)
        {
            fprintf(magnitudeFile, "%u, %d\n", i, ppg_buffer[i]);
            fflush(magnitudeFile);
        }
    }
#endif

    //Step 1: Apply band pass filter
    bandpass_buffer(ppg_buffer, filtered_buff);

    //Step 2: Remove the mean
    remove_mean(filtered_buff);

    //Step 3: Calculate autocorrelation
    autocorr(filtered_buff, autocorr_buff);


    //Next step: Simply find the peak of the autocorrelation signal
    uint8_t peak_ind = 0;
    float interval = 0;
    uint8_t i;
    for (i = FIRST_AUTOCORR_PEAK_LAG; i < NUM_AUTOCORR_LAGS; i++)
    {

        if ((autocorr_buff[i] > autocorr_buff[i - 1]) && (autocorr_buff[i] > autocorr_buff[i + 1])){
            peak_ind = i;
            break;
        }
        
    }

    //Skip this for now
    /*
    // Step 7: Check the conditions to see if it was a real peak or not, and if so, keep the peak as valid
    uint8_t neg_slope_count = 0;
    int64_t delta_amplitude_right = 0;
    uint8_t pos_slope_count = 0;
    int64_t delta_amplitude_left = 0;
    get_autocorr_peak_stats(autocorr_buff, &neg_slope_count, &delta_amplitude_right, &pos_slope_count, &delta_amplitude_left, peak_ind);
    if ((pos_slope_count > AUTOCORR_MIN_HALF_LEN) && (neg_slope_count > AUTOCORR_MIN_HALF_LEN) && (delta_amplitude_right > AUTOCORR_DELTA_AMPLITUDE_THRESH) && (delta_amplitude_left > AUTOCORR_DELTA_AMPLITUDE_THRESH))
    {
        // The period is peak_ind/sampling_rate seconds. That corresponds to a frequency of 1/period
        // Then the HR is computed doing 60/period
        interval = ((float)(peak_ind) / (float)(SAMPLING_RATE));
    }
    else
    {
        // not a valid autocorrelation peak
        interval = 0;
    }
        */


    if (peak_ind != 0) {
        interval = ((float)(peak_ind) / (float)(SAMPLING_RATE));
    } else {
        interval = 0;
    }
 
    return interval;

}

//Final function called in the main script:
int autocorrelation_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz)
{
    //Saving the samples in the buffer
    autocorr_buffer[autocorr_buffer_index] = ppg;
    autocorr_buffer_index++;
    float temp = 0;

    if (autocorr_buffer_index >= NUM_TUPLES) //When the window length is reached
    {
        temp = autcorr_count_beats(autocorr_buffer);
        autocorr_buffer_index = 0;
    }

    if (temp != 0) //It means the autocorrelation peak has been found
    {
        autocorrelation_HR = (int)(600/temp); // Temp is the interval betweeen two peaks, here: (1/T)*60*10 to obtain bpm*10
        
        //MEDIAN FILTER on the result
        //Saving the results in a buffer to implement a moving average or median
        autocorr_results[autocorr_results_index] = autocorrelation_HR; 
        autocorr_results_index++;
        if (autocorr_results_index >= HRM_HIST_LEN)
        {
            autocorr_results_index = 0;
        }
        //Median
        //First: Sorting the buffer
        bool busy;
        do {
            busy = false;
            for (int i=0;i<HRM_HIST_LEN-1;i++) {
                if (autocorr_results[i] > autocorr_results[i+1]) {
                    int te = autocorr_results[i];
                    autocorr_results[i] = autocorr_results[i+1];
                    autocorr_results[i+1] = te;
                    busy = true;
                }
            }
        } while (busy);
        // calculate HRM from middle values
        int16_t min = (int16_t)((HRM_HIST_LEN - HRM_MEDIAN_LEN)/2); // 1
        int16_t max = (int16_t)((HRM_HIST_LEN + HRM_MEDIAN_LEN)/2); // 3
        int16_t n = 0;
        int sumBPM = 0;
        for (int16_t i=min;i<max;i++) {
            if (autocorr_results[i]==0) continue;
            sumBPM += autocorr_results[i];
            n++;
        }
        if (n) {
            autocorrelation_HR = (int)(sumBPM/n);
        } 
    }
    printf("HR: %d\n", autocorrelation_HR);
    return autocorrelation_HR; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Second repiclica applying the peak verification and the derivative signal

/*

#include "autocorrelationHeartRate.h"
#include <stdint.h>
#include <stdio.h>  
#include <string.h>
#include <stdbool.h>

#include "../../utils/bandpass_filter/bandpass_filter.h"
#include "../../types.h"

#define SAMPLING_RATE 25  // 25 hz sampling rate
#define NUM_TUPLES 125    // 5 seconds worth of data
#define WINDOW_LENGTH (NUM_TUPLES / SAMPLING_RATE) // window length in seconds: 5

static int autocorrelation_HR; 
static int autocorrelation_HR_temp = 0;

//Da togliere:
//static ppg_t autocorr_buffer[NUM_TUPLES]; //Buffer where I put the initial signal
//static long autocorr_buffer_index = 0;


#define NUM_AUTOCORR_LAGS 51      // number of lags to calculate for autocorrelation. At a minimum heart rate of 30 bpm -> 0.5 beat/s -> 2 / sampling_period(0.04s) = 2s * 25Hz = 50 lags 
#define FIRST_AUTOCORR_PEAK_LAG 8 // corresponds to the first feasible autocorrelation lag -> at a max heart rate of 3 beats /s -> 0.333s / sampling_period(0.04s) = 8.25 -> 8 lags

#define DERIV_FILT_LEN 7                    // length of derivative filter
#define AUTOCORR_DELTA_AMPLITUDE_THRESH 1e4 // this is the min delta between peak and trough of autocorrelation peak
#define AUTOCORR_MIN_HALF_LEN 2             // this is the min number of points the autocorrelation peak should be on either side of the peak

//DUMP FILE: used to print the signal in .csv files
//#define DUMP_FILE 
#ifdef DUMP_FILE
static int autocorr_passes = 0; // counter of how many times the autocorr has been called
#define DUMP_REMOVED_MEAN_FILE_NAME "removed_mean.csv"
#define DUMP_AUTOCORRELATION_FILE_NAME "autocorrelation"
#define DUMP_MAGNITUDE_FILE_NAME "magnitude.csv"
#define DUMP_FILTERED_FILE_NAME "filtered.csv"
static FILE *autocorrelationFile;
static FILE *magnitudeFile;
static FILE *filteredFile;
static FILE *removedMeanFile;
static FILE *autocorrelationFile;
#endif

static ppg_t autocorr_buffer[NUM_TUPLES];
static long autocorr_buffer_index = 0;
static ppg_t filtered_buff[NUM_TUPLES] = {0};    // low pass filtered data
static int64_t autocorr_buff[NUM_AUTOCORR_LAGS] = {0}; // autocorrelation results
static int64_t deriv[NUM_AUTOCORR_LAGS] = {0};         // derivative
//Buffer to implement median filter on the estimated HR
#define HRM_HIST_LEN 4
#define HRM_MEDIAN_LEN 2
static int autocorr_results[HRM_HIST_LEN] = {0};
static int autocorr_results_index = 0;

//BandPass filter
static BPFilter bpFilter;

void autocorrelation_heartrate_init()
{
    BPFilter_init(&bpFilter);
    //Putting everything to zero
    for (int i = 0; i < NUM_AUTOCORR_LAGS; i++)
    {
        autocorr_buff[i] = 0;
        deriv[i] = 0;
    }
    autocorrelation_HR = 0;
    for (int i = 0; i < NUM_TUPLES; i++)
    {
        autocorr_buffer[i] = 0;;
    }
    autocorr_buffer_index = 0;
    for (int i = 0; i < HRM_HIST_LEN; i++)
    {
        autocorr_results[i] = 0;;
    }
    autocorr_results_index = 0;


#ifdef DUMP_FILE
    magnitudeFile = fopen(DUMP_MAGNITUDE_FILE_NAME, "w");
    filteredFile = fopen(DUMP_FILTERED_FILE_NAME, "w");
    removedMeanFile = fopen(DUMP_REMOVED_MEAN_FILE_NAME, "w");

    autocorr_passes = 0;
#endif

}   

// Sends the buffer to a band pass filter
static void bandpass_buffer(ppg_t *input_buffer, ppg_t *output_buffer)
{
    for (int i = 0; i < NUM_TUPLES; i++)
    {
        BPFilter_put(&bpFilter, input_buffer[i]);
        output_buffer[i] = BPFilter_get(&bpFilter);

#ifdef DUMP_FILE
        if (filteredFile)
        {
            fprintf(filteredFile, "%d\n", output_buffer[i]);
            fflush(filteredFile);
        }
#endif
    }
}

// Remove mean from filtered data
static void remove_mean(ppg_t *buffer)
{
    int64_t sum_temp = 0; //The sum_temp is added to be sure to convert correctly between different data type (probably it is not necessary)
    ppg_t sum = 0;
    uint16_t i;
    for (i = 0; i < NUM_TUPLES; i++)
    {
        sum_temp += (int64_t)buffer[i];
    }
    sum_temp = (int64_t)((float)(sum_temp) / (float)(NUM_TUPLES));
    sum = (ppg_t)sum_temp;
    for (i = 0; i < NUM_TUPLES; i++)
    {
        buffer[i] -= sum;

#ifdef DUMP_FILE
        if (removedMeanFile)
        {
            fprintf(removedMeanFile, "%d\n", buffer[i]);
        }
#endif
    }

#ifdef DUMP_FILE
    if (removedMeanFile)
        fflush(removedMeanFile);
#endif
}

// Autocorrelation function: it computes the autocorrelatio for lags between 0 and NUM_AUTOCORR_LAGS
static void autocorr(ppg_t *buffer, int64_t *autocorr_buff)
{

#ifdef DUMP_FILE
    autocorr_passes++;
    char autocorrFileName[100] = DUMP_AUTOCORRELATION_FILE_NAME;
    char idxstr[5];
    sprintf(idxstr, "%d", autocorr_passes);
    strcat(autocorrFileName, idxstr);
    strcat(autocorrFileName, ".csv");
    autocorrelationFile = fopen(autocorrFileName, "w+");
#endif


    uint8_t lag;
    uint16_t i;
    int64_t temp_ac;
    for (lag = 0; lag < NUM_AUTOCORR_LAGS; lag++)
    {
        temp_ac = 0;
        for (i = 0; i < NUM_TUPLES - lag; i++)
        {
            temp_ac += (int64_t)buffer[i] * (int64_t)buffer[i + lag];
        }
        autocorr_buff[lag] = temp_ac;

#ifdef DUMP_FILE
        if (autocorrelationFile)
        {
            fprintf(autocorrelationFile, "%u, %lld\n", lag, autocorr_buff[lag]);
        }
#endif
    }
#ifdef DUMP_FILE
    if (autocorrelationFile)
    {
        fflush(autocorrelationFile);
        fclose(autocorrelationFile);
    }
#endif
}

// Derivative calculation
static void derivative(int64_t *autocorr_buff, int64_t *deriv)
{
    uint8_t n = 0;
    uint8_t i = 0;
    int64_t temp_deriv = 0;

    //If I use the derivative filters with coefficients (not used in this mplementation)
    //     for (n = 0; n < NUM_AUTOCORR_LAGS; n++)
    //     {
    //         temp_deriv = 0;
    //         for (i = 0; i < DERIV_FILT_LEN; i++)
    //         {
    //             if (n - i >= 0)
    //             {
    //                 temp_deriv += deriv_coeffs[i] * autocorr_buff[n - i];
    //             }
    //         }
    //         deriv[n] = temp_deriv;


    for (n = 0; n < NUM_AUTOCORR_LAGS; n++)
    {
        if (n > 0)
        {
            deriv[n] = autocorr_buff[n] - autocorr_buff[n - 1];
        }
        else
            deriv[0] = 0;

    }   
}

//function to get the precise peak: if the found peak is not actually the peak -> Moving to left or right to find the actual peak
static uint8_t get_precise_peakind(int64_t *autocorr_buff, uint8_t peak_ind)
{
    uint8_t loop_limit = 0;
    if ((autocorr_buff[peak_ind] > autocorr_buff[peak_ind - 1]) && (autocorr_buff[peak_ind] > autocorr_buff[peak_ind + 1]))
    {
        // peak_ind is perfectly set at the peak. nothing to do
    }
    else if ((autocorr_buff[peak_ind] > autocorr_buff[peak_ind + 1]) && (autocorr_buff[peak_ind] < autocorr_buff[peak_ind - 1]))
    {
        // peak is to the left. keep moving in that direction
        loop_limit = 0;
        while ((autocorr_buff[peak_ind] > autocorr_buff[peak_ind + 1]) && (autocorr_buff[peak_ind] < autocorr_buff[peak_ind - 1]) && (loop_limit < 10))
        {
            peak_ind = peak_ind - 1;
            loop_limit++;
        }
    }
    else
    {
        // peak is to the right. keep moving in that direction
        loop_limit = 0;
        while ((autocorr_buff[peak_ind] > autocorr_buff[peak_ind - 1]) && (autocorr_buff[peak_ind] < autocorr_buff[peak_ind + 1]) && (loop_limit < 10))
        {
            peak_ind = peak_ind + 1;
            loop_limit++;
        }
    }
    return peak_ind;
}

// take a look at the original autocorrelation signal at index i and see if
// it's a real peak or if it's just a fake "noisy" peak corresponding to
// non-walking. Basically just count the number of points of the
// autocorrelation peak to the right and left of the peak. this function gets
// the number of points to the right and left of the peak, as well as the delta amplitude
static void get_autocorr_peak_stats(int64_t *autocorr_buff, uint8_t *neg_slope_count, int64_t *delta_amplitude_right, uint8_t *pos_slope_count, int64_t *delta_amplitude_left, uint8_t peak_ind)
{

    // first look to the right of the peak. walk forward until the slope begins decreasing (increasing in my opinion)
    uint8_t neg_slope_ind = peak_ind;
    uint16_t loop_limit = NUM_AUTOCORR_LAGS - 1;
    while ((autocorr_buff[neg_slope_ind + 1] - autocorr_buff[neg_slope_ind] < 0) && (neg_slope_ind < loop_limit))
    {
        *neg_slope_count = *neg_slope_count + 1;
        neg_slope_ind = neg_slope_ind + 1;
    }

    // get the delta amplitude between peak and right trough (negative peak)
    *delta_amplitude_right = autocorr_buff[peak_ind] - autocorr_buff[neg_slope_ind];

    // next look to the left of the peak. walk backward until the slope begins increasing (decreasing in my opinion)
    uint8_t pos_slope_ind = peak_ind;
    loop_limit = 0;
    while ((autocorr_buff[pos_slope_ind] - autocorr_buff[pos_slope_ind - 1] > 0) && (pos_slope_ind > loop_limit))
    {
        *pos_slope_count = *pos_slope_count + 1;
        pos_slope_ind = pos_slope_ind - 1;
    }

    // get the delta amplitude between the peak and the left trough
    *delta_amplitude_left = autocorr_buff[peak_ind] - autocorr_buff[pos_slope_ind];
}

// Function that returns the period in seconds of the signal (interval) -> it corresponds to the lag of the second peak of autocorrelation divided by f_sample
static float autcorr_count_beats(ppg_t *ppg_buffer)
{

#ifdef DUMP_FILE
    for (int i = 0; i < NUM_TUPLES; i++)
    {
        if (magnitudeFile)
        {
            fprintf(magnitudeFile, "%u, %d\n", i, ppg_buffer[i]);
            fflush(magnitudeFile);
        }
    }
#endif

    //Step 1: Apply band pass filter
    bandpass_buffer(ppg_buffer, filtered_buff);

    //Step 2: Remove the mean
    remove_mean(filtered_buff);

    //Step 3: Calculate autocorrelation
    autocorr(filtered_buff, autocorr_buff);

    // Step 4: Calculate derivative
    derivative(autocorr_buff, deriv);

    // Step 5: find peak
    // look for first zero crossing where derivative goes from positive to negative. That corresponds to the first positive peak in the autocorrelation.
    // Look at two samples instead of just one to maybe reduce the chances of getting tricked by noise
    uint8_t peak_ind = 0;
    float interval = 0;
    uint8_t i;
    for (i = FIRST_AUTOCORR_PEAK_LAG; i < NUM_AUTOCORR_LAGS; i++)
    {
        if ((deriv[i] < 0) && (deriv[i - 1] < 0) && (deriv[i - 2] > 0) && (deriv[i - 3] > 0))
        {
            peak_ind = i - 1;
            break;
        } 
    }

    // Step 6: Hone in on the exact peak index
    peak_ind = get_precise_peakind(autocorr_buff, peak_ind);
    // printf("peak ind: %i\n", peak_ind);


    // Step 7: Check the conditions to see if it was a real peak or not, and if so, keep the peak as valid
    uint8_t neg_slope_count = 0;
    int64_t delta_amplitude_right = 0;
    uint8_t pos_slope_count = 0;
    int64_t delta_amplitude_left = 0;
    get_autocorr_peak_stats(autocorr_buff, &neg_slope_count, &delta_amplitude_right, &pos_slope_count, &delta_amplitude_left, peak_ind);
    if ((pos_slope_count > AUTOCORR_MIN_HALF_LEN) && (neg_slope_count > AUTOCORR_MIN_HALF_LEN) && (delta_amplitude_right > AUTOCORR_DELTA_AMPLITUDE_THRESH) && (delta_amplitude_left > AUTOCORR_DELTA_AMPLITUDE_THRESH))
    {
        // The period is peak_ind/sampling_rate seconds. That corresponds to a frequency of 1/period
        // Then the HR is computed doing 60/period
        interval = ((float)(peak_ind) / (float)(SAMPLING_RATE));
    }
    else
    {
        // not a valid autocorrelation peak
        interval = 0;
    }

 
    return interval;
    
}

//Final function called in the main script:
int autocorrelation_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz)
{
    //Saving the samples in the buffer
    autocorr_buffer[autocorr_buffer_index] = ppg;
    autocorr_buffer_index++;
    float temp = 0;

    if (autocorr_buffer_index >= NUM_TUPLES) //When the window length is reached
    {
        temp = autcorr_count_beats(autocorr_buffer);
        autocorr_buffer_index = 0;
    }

    if (temp != 0) //It means the autocorrelation peak has been found
    {
        autocorrelation_HR = (int)(600/temp); // Temp is the interval betweeen two peaks, here: (1/T)*60*10 to obtain bpm*10
        
        //MEDIAN FILTER on the result
        //Saving the results in a buffer to implement a moving average or median
        autocorr_results[autocorr_results_index] = autocorrelation_HR; 
        autocorr_results_index++;
        if (autocorr_results_index >= HRM_HIST_LEN)
        {
            autocorr_results_index = 0;
        }
        //Median
        //First: Sorting the buffer
        bool busy;
        do {
            busy = false;
            for (int i=0;i<HRM_HIST_LEN-1;i++) {
                if (autocorr_results[i] > autocorr_results[i+1]) {
                    int te = autocorr_results[i];
                    autocorr_results[i] = autocorr_results[i+1];
                    autocorr_results[i+1] = te;
                    busy = true;
                }
            }
        } while (busy);
        // calculate HRM from middle values
        int16_t min = (int16_t)((HRM_HIST_LEN - HRM_MEDIAN_LEN)/2); // 1
        int16_t max = (int16_t)((HRM_HIST_LEN + HRM_MEDIAN_LEN)/2); // 3
        int16_t n = 0;
        int sumBPM = 0;
        for (int16_t i=min;i<max;i++) {
            if (autocorr_results[i]==0) continue;
            sumBPM += autocorr_results[i];
            n++;
        }
        if (n) {
            autocorrelation_HR = (int)(sumBPM/n);
        } 
    } 
    printf("HR: %d\n", autocorrelation_HR);
    return autocorrelation_HR; 
}

*/