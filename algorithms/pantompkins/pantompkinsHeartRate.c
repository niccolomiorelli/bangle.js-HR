/* ----------------------------------------------------
 * PANTOMPKINS ALGORITHM
 * ----------------------------------------------------
 * ------------------------------------------------------------------------------*
 * File: panTompkins.c                                                           *
 *       ANSI-C implementation of Pan-Tompkins real-time QRS detection algorithm *
 * Author: Rafael de Moura Moreira <rafaelmmoreira@gmail.com>                    *
 * License: MIT License                                                          *
 * ------------------------------------------------------------------------------*
 * ---------------------------------- HISTORY ---------------------------------- *
 *    date   |    author    |                     description                    *
 * ----------| -------------| ---------------------------------------------------*
 * 2019/04/11| Rafael M. M. | - Fixed moving-window integral.                    *
 *           |              | - Fixed how to find the correct sample with the    *
 *           |              | last QRS.                                          *
 *           |              | - Replaced constant value in code by its #define.  *
 *           |              | - Added some casting on comparisons to get rid of  *
 *           |              | compiler warnings.                                 *
 * 2019/04/15| Rafael M. M. | - Removed delay added to the output by the filters.*
 *           |              | - Fixed multiple detection of the same peak.       *
 * 2019/04/16| Rafael M. M. | - Added output buffer to correctly output a peak   *
 *           |              | found by back searching using the 2nd thresholds.  *
 * 2019/04/23| Rafael M. M. | - Improved comparison of slopes.                   *
 *           |              | - Fixed formula to obtain the correct sample from  *
 *           |              | the buffer on the back search.                     *
 * ------------------------------------------------------------------------------*
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2018 Rafael de Moura Moreira                                    *
 *                                                                               *
 * Permission is hereby granted, free of charge, to any person obtaining a copy  *
 * of this software and associated documentation files (the "Software"), to deal *
 * in the Software without restriction, including without limitation the rights  *
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     *
 * copies of the Software, and to permit persons to whom the Software is         *
 * furnished to do so, subject to the following conditions:                      *
 *                                                                               *
 * The above copyright notice and this permission notice shall be included in all*
 * copies or substantial portions of the Software.                               *
 *                                                                               *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, *
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE *
 * SOFTWARE.                                                                     *
 *-------------------------------------------------------------------------------*
 * Description                                                                   *
 *                                                                               *
 * The main goal of this implementation is to be easy to port to different opera-*
 * ting systems, as well as different processors and microcontrollers, including *
 * embedded systems. It can work both online or offline, depending on whether all*
 * the samples are available or not - it can be adjusted on the input function.  *
 *                                                                               *
 * ----------------------------------------------------------------------------- *
 *                                                                               *
 * The algoritm has been modified and adapted to address PPG signals. Here the   *
 * modified version is described.                                                *
 *                                                                               *
 * The main function, pantompkings_heartrate(), receives as input the new sample *
 * and stores it in a buffer. Then it runs through a chain of filters: DC block, *
 * bandpass filter (0.5 - 4 Hz), derivative filter, squaring and windowed inte-  *
 * gratior. 
 * 
 * Then, signal peaks in the integrator signal are identified, and to be consi-  *
 * dered as an actual peak they must be above a certain threshold (threshold_i1) *
 * and the corresponding value in the filtered ppg signal must be above a second *                                                                        
 * threshold (threshold_f1). The thresholds on th integrated signal is updated   *
 * following the classical Pan-Tompkins algorithm in this way:                   *
 * [threshold_i1 = npk_i + 0.25*(spk_i - npk_i)].  Every peak, whether it is a   *
 * an actual peak or noise contributes to update this threshold. While the thre- *
 * shold on the filtered signal is an exponential moving average of the filtered *
 * signal.                                                                       *                     
 * Additionally, the time-restraints normally used in Pan-Tompkins are kept.     *
 * A hard 200ms restrain (a new peak within 200ms is considered noise) and a soft*
 * 360ms restrain (the peak's squared slope must also be very high to be consi-  *
 * dered as a real peak).                                                        *                
 *                                                                               *
 * Two buffers keep 8 RR-intervals to calculate RR-averages: one of them keeps   *
 * the last 8 RR-intervals, while the other keeps only the RR-intervals that res-*
 * pect certain restrictions. If both averages are equal, the heart pace is con- *
 * sidered normal. If the heart rate isn't normal, the thresholds change to make *
 * it easier to detect possible weaker peaks. If no peak is detected for a long  *
 * period of time, the thresholds also change and the last discarded peak candi- *
 * date is reconsidered.                                                         *
 *                                                                               *
 * A median Filter is implemented on the output HR, to discard outliers and      *
 * smooth the output.                                                            *
 *                                                                               *
 * NOTE: Possible further implementation:                                        *
 * A back search can be implemented to look for missed peaks.                    *
 * When a HR value is abnormal (outlier due to noise missdetection of peaks),    *
 * the algorithm can give rravvg2 (avarge of good RR-intervas), but for now this *
 * is not implemented because there could be problems during the starting that   *
 * can affect the whole algorithm.                                               *  
 *-------------------------------------------------------------------------------*
*/

#include <stdio.h>
#include <stdbool.h>


#include "../../utils/bandpass_filter/bandpass_filter.h"
#include "../../types.h"

#include "pantompkinsHeartRate.h"

#define WINDOWSIZE 5   // Integrator window size, in samples. Around 150 ms is suggested
						
#define NOSAMPLE -32000 // An indicator that there are no more samples to read. Use an impossible value for a sample.
#define FS 25          // Sampling frequency.
#define BUFFSIZE 40    // The size of the buffers (in samples). Still to decide                         

#define DELAY 17 //Delay in samples introduced by the filter 

//For file dumping
// #define DUMP_FILE
#ifdef DUMP_FILE
#include <stdio.h>
#define DUMP_RAW_FILE_NAME "raw.csv"
#define DUMP_DCREMOVE_FILE_NAME "dcremove.csv"
#define DUMP_FILTERED_FILE_NAME "filtered.csv"
#define DUMP_DERIVATIVE_FILE_NAME "derivative.csv"
#define DUMP_SQUARED_FILE_NAME "squared.csv"
#define DUMP_INTEGRATED_FILE_NAME "integrated.csv"
#define DUMP_THRESHOLDS_FILE_NAME "thresholds.csv"
#define DUMP_QRS_FILE_NAME "qrs.csv"

static FILE *rawFile;
static FILE *dcremoveFile;
static FILE *filteredFile;
static FILE *derivativeFile;
static FILE *squaredFile;
static FILE *integratedFile;
static FILE *thresholdsFile;
static FILE *qrsFile;
#endif

//Defining the variables
// The signal array is where the most recent samples are kept. The other arrays are the outputs of each
// filtering module: DC Block, low pass, high pass, integral etc.
// The output is a buffer where we can change a previous result (using a back search) before outputting.
static ppg_t signal[BUFFSIZE], dcblock[BUFFSIZE], bandpass[BUFFSIZE], derivative[BUFFSIZE], squared[BUFFSIZE], integral[BUFFSIZE], outputSignal[BUFFSIZE];

// rr1 holds the last 8 RR intervals. rr2 holds the last 8 RR intervals between rrlow and rrhigh.
// rravg1 is the rr1 average, rr2 is the rravg2. rrlow = 0.92*rravg2, rrhigh = 1.08*rravg2 and rrmiss = 1.16*rravg2.
// rrlow is the lowest RR-interval considered normal for the current heart beat, while rrhigh is the highest.
// rrmiss is the longest that it would be expected until a new QRS is detected. If none is detected for such
// a long interval, the thresholds must be adjusted.
static int rr1[8], rr2[8], rravg1, rravg2, rrlow = 0, rrhigh = 0, rrmiss = 0;

// i and j are iterators for loops.
// sample counts how many samples have been read so far.
// lastQRS stores which was the last sample read when the last R sample was triggered.
// lastSlope stores the value of the squared slope when the last R sample was triggered.
// currentSlope helps calculate the max. square slope for the present sample.
// These are all long unsigned int so that very long signals can be read without messing the count.
static long unsigned int i, j, sample = 0, lastQRS = 0, lastSlope = 0, currentSlope = 0;

// This variable is used as an index to work with the signal buffers. If the buffers still aren't
// completely filled, it shows the last filled position. Once the buffers are full, it'll always
// show the last position, and new samples will make the buffers shift, discarding the oldest
// sample and storing the newest one on the last position.
static int current;

// There are the variables from the original Pan-Tompkins algorithm.
// The ones ending in _i correspond to values from the integrator.
// The ones ending in _f correspond to values from the DC-block/low-pass/high-pass filtered signal.
// The peak variables are peak candidates: signal values above the thresholds.
// The threshold 1 variables are the threshold variables. If a signal sample is higher than this threshold, it's a peak.
// The threshold 2 variables are half the threshold 1 ones. They're used for a back search when no peak is detected for too long.
// The spk and npk variables are, respectively, running estimates of signal and noise peaks.
static ppg_t peak_i = 0, peak_f = 0, threshold_i1 = 0, threshold_i2 = 0, threshold_f1 = 0, threshold_f2 = 0, spk_i = 0, spk_f = 0, npk_i = 0, npk_f = 0;

// qrs tells whether there was a detection or not.
// regular tells whether the heart pace is regular or not.
// prevRegular tells whether the heart beat was regular before the newest RR-interval was calculated.
static bool qrs, regular = true, prevRegular;

//Band Pass Filter
static BPFilter bpFilter;

//Added variables
static bool peak = false;
static num_beats = 0;
int HR;

//For the median filter
//Buffer per implementare media mobile o mediana sui risultati
#define HRM_HIST_LEN 8
#define HRM_MEDIAN_LEN 4
static int pt_results[HRM_HIST_LEN] = {0};
static int pt_results_index = 0;

void pantompkins_heartrate_init(){

    //Initializations:
    rravg1, rravg2, rrlow = 0;
    rrhigh = 0, rrmiss = 0;
    i, j, sample = 0, lastQRS = 0;
    lastSlope = 0;
    currentSlope = 0;
    peak_i = 0, peak_f = 0, threshold_i1 = 0, threshold_i2 = 0, threshold_f1 = 0, threshold_f2 = 0, spk_i = 0, spk_f = 0, npk_i = 0, npk_f = 0;
    sample = 0;
    regular = true;
    for (i=0; i < BUFFSIZE; i++)
    {
        signal[i] = 0;
        dcblock[i] = 0;
        bandpass[i] = 0;
        derivative[i] = 0;
        squared[i] = 0;
        integral[i] = 0;
        outputSignal[i] = 0;
    }

    // Initializing the RR averages
	for (i = 0; i < 8; i++)
    {
        rr1[i] = 0;
        rr2[i] = 0;
    }

    //Bandpass filter
    BPFilter_init(&bpFilter);
    

    peak = false;
    lastQRS = 0;
    num_beats = 0;

    for (int i = 0; i < HRM_HIST_LEN; i++)
    {
        pt_results[i] = 0;;
    }
    pt_results_index = 0;

#ifdef DUMP_FILE
    rawFile = fopen(DUMP_RAW_FILE_NAME, "w+");
    dcremoveFile = fopen(DUMP_DCREMOVE_FILE_NAME, "w+");
    filteredFile = fopen(DUMP_FILTERED_FILE_NAME, "w+");
    derivativeFile = fopen(DUMP_DERIVATIVE_FILE_NAME, "w+");
    squaredFile = fopen(DUMP_SQUARED_FILE_NAME, "w+");
    integratedFile = fopen(DUMP_INTEGRATED_FILE_NAME, "w+");
    thresholdsFile = fopen(DUMP_THRESHOLDS_FILE_NAME, "w+");
    qrsFile = fopen(DUMP_QRS_FILE_NAME, "w+");
#endif


}

int pantompkins_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz){
    // Test if the buffers are full.
    // If they are, shift them, discarding the oldest sample and adding the new one at the end.
    // Else, just put the newest sample in the next free position.
    // Update 'current' so that the program knows where's the newest sample.
    if (sample >= BUFFSIZE)
    {
        for (i = 0; i < BUFFSIZE - 1; i++)
        {
            signal[i] = signal[i+1];
            dcblock[i] = dcblock[i+1];
            bandpass[i] = bandpass[i+1];
            derivative[i] = derivative[i+1];
            squared[i] = squared[i+1];
            integral[i] = integral[i+1];
            outputSignal[i] = outputSignal[i+1];
        }
        current = BUFFSIZE - 1;
    }
    else
    {
        current = sample;
    }
    signal[current] = ppg; //Adding the new sample to the signal buffer
    sample++; // Update sample counter

#ifdef DUMP_FILE
        if (rawFile)
        {
            if (!fprintf(rawFile, "%d,\n", signal[current]))
                puts("error writing file");
            fflush(rawFile);
        }
#endif

    //STEP 1: DC BLOCKING
    // DC Block filter
    // This was not proposed on the original paper.
    // It is not necessary and can be removed if your sensor or database has no DC noise.
    if (current >= 1)
        dcblock[current] = signal[current] - signal[current-1] + 0.995*dcblock[current-1]; //NOTA: current punta all'ultio elemento del buffer, quindi l'ultimo sample
    else
        dcblock[current] = 0;

#ifdef DUMP_FILE
    if (dcremoveFile)
    {
        if (!fprintf(dcremoveFile, "%d,\n", dcblock[current]))
            puts("error writing file");
        fflush(dcremoveFile);
    }
#endif

    //STEP 2: BANDPASS FILTER
    // Applying the filter
    BPFilter_put(&bpFilter,dcblock[current]);
    bandpass[current] = BPFilter_get(&bpFilter);

#ifdef DUMP_FILE
    if (filteredFile)
    {
        if (!fprintf(filteredFile, "%d,\n", bandpass[current]))
            puts("error writing file");
        fflush(filteredFile);
    }
#endif

    //STEP 3: DERIVATIVE
    derivative[current] = bandpass[current];
    if (current > 0)
        derivative[current] -= bandpass[current-1];

#ifdef DUMP_FILE
    if (derivativeFile)
    {
        if (!fprintf(derivativeFile, "%d,\n", derivative[current]))
            puts("error writing file");
        fflush(derivativeFile);
    }
#endif

    //STEP 4: SQUARING
    squared[current] = derivative[current]*derivative[current];

#ifdef DUMP_FILE
    if (squaredFile)
    {
        if (!fprintf(squaredFile, "%d,\n", squared[current]))
            puts("error writing file");
        fflush(squaredFile);
    }
#endif


    //STEP 5: INTEGRATION
    integral[current] = 0;
    for (i = 0; i < WINDOWSIZE; i++)
    {
        if (current >= (ppg_t)i)
            integral[current] += squared[current - i];
        else
            break;
    }
    integral[current] /= (ppg_t)i;

#ifdef DUMP_FILE
    if (integratedFile)
    {
        if (!fprintf(integratedFile, "%d,\n", integral[current]))
            puts("error writing file");
        fflush(integratedFile);
    }
#endif



    qrs = false;

    //Algorithm form now on:
 
    //Updating the threshold for the bandpass filter -> Simply a moving average of the signal with a quite reactive factor (alpha = 0.1)
    //It is adopted because the filtered PPG can have a fluctauting baseline, so I need a reactive threshold
    threshold_f1 = 0.1*bandpass[current] + 0.9*threshold_f1; 

    //Looking for peaks in the integrated signal, checking that those point are above threshold in the bandpass filtered signal
    if(sample>=3 && integral[current-1]>=integral[current-2] && integral[current-1]>=integral[current] && bandpass[current-1]>= threshold_f1){
        peak = true;
        peak_i = integral[current-1];
    }else peak = false;


    //If peaks in the integral signal are also above threshold, they're probably signal peaks.
    if ((integral[current-1] >= threshold_i1 && peak)) 
    {
        // There's a 200ms latency. If the new peak respects this condition, we can keep testing.
        if (sample > lastQRS + FS*0.3)
        {
            // If it respects the 200ms latency, but it doesn't respect the 360ms latency, we check the slope.
            if (sample <= lastQRS + (long unsigned int)(0.36*FS))
            {
                // The squared slope is "M" shaped. So we have to check nearby samples to make sure we're really looking
                // at its peak value, rather than a low one.
                currentSlope = 0; //Squared signal value
                for (j = current - 10; j <= current; j++){
                    if (squared[j] > currentSlope)
                        currentSlope = squared[j];
                }

                if (currentSlope <= (ppg_t)(lastSlope/2))
                {
                    qrs = false;
                } else
                {
                    spk_i = 0.125*peak_i + 0.875*spk_i;
                    threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
                    threshold_i2 = 0.5*threshold_i1;

                    spk_f = 0.125*peak_f + 0.875*spk_f;

                    lastSlope = currentSlope;
                    qrs = true;
                }
            }
            // If it was above both thresholds and respects both latency periods, it certainly is a R peak.
            else
            {
                currentSlope = 0;
                for (j = current - 10; j <= current; j++)
                    if (squared[j] > currentSlope)
                        currentSlope = squared[j];

                spk_i = 0.125*peak_i + 0.875*spk_i;
                threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
                threshold_i2 = 0.5*threshold_i1;

                spk_f = 0.125*peak_f + 0.875*spk_f;

                lastSlope = currentSlope;
                qrs = true;
            }
        }
        // If the new peak doesn't respect the 200ms latency, it's noise. Update thresholds and move on to the next sample.
        else
        {
            npk_i = 0.125*peak_i + 0.875*npk_i;
            threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
            threshold_i2 = 0.5*threshold_i1;
            peak_f = bandpass[current];
            qrs = false;
        }
        
    } else if (integral[current-1]<=threshold_i1 && peak){ //If we have a peak but it's lower than the threshold
        npk_i = 0.125*peak_i + 0.875*npk_i;
        threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
        threshold_i2 = 0.5*threshold_i1;
        peak_f = bandpass[current];
        qrs = false;
    }

    
#ifdef DUMP_FILE
    if (thresholdsFile)
    {
        if (!fprintf(thresholdsFile, "%d,%d, %d, %d \n", threshold_f1,threshold_i1,npk_i,spk_i))
            puts("error writing file");
        fflush(thresholdsFile);
    }
#endif

    //--------------------------------------------------- Second part of the algorithm

    // If a R-peak was detected, the RR-averages must be updated.
    if (qrs)
    {
        num_beats++;
        // Add the newest RR-interval to the buffer and get the new average.
        rravg1 = 0;
        for (i = 0; i < 7; i++)
        {
            rr1[i] = rr1[i+1];
            rravg1 += rr1[i];
        }
        rr1[7] = sample - lastQRS;
        lastQRS = sample;
        rravg1 += rr1[7];
        rravg1 *= 0.125;

        // If the newly-discovered RR-average is normal, add it to the "normal" buffer and get the new "normal" average.
        // Update the "normal" beat parameters.
        if (num_beats <= 8) {
            rr2[num_beats] = rr1[7];
            num_beats++;
            rravg2 = 0;
            for (i = 0; i < num_beats; i++)
                rravg2 += rr2[i];
            rravg2 /= num_beats;
            rrlow = 0.7 * rravg2;
            rrhigh = 1.3 * rravg2;
            rrmiss = 1.75 * rravg2;
            //Logic in case of problems: Considering physiological values
            if (rrlow < 7) rrlow = 7; //285ms intra beat interval -> 210 bpm
            if (rrhigh > 50) rrhigh = 50; //2s intra beat interval -> 30 bpm
            if (rrmiss > 90) rrmiss = 90; //it means no beat for 3.6s (17bpm allowd, clearly out of range)
        } else if ( (rr1[7] >= rrlow) && (rr1[7] <= rrhigh) )
        {
            rravg2 = 0;
            for (i = 0; i < 7; i++)
            {
                rr2[i] = rr2[i+1];
                rravg2 += rr2[i];
            }
            rr2[7] = rr1[7];
            rravg2 += rr2[7];
            rravg2 *= 0.125;
            rrlow = 0.7*rravg2;
            rrhigh = 1.3*rravg2;
            rrmiss = 1.75*rravg2;
            //Logic in case of problems: Considering physiological values
            if (rrlow < 7) rrlow = 7; //285ms intra beat interval -> 210 bpm
            if (rrhigh > 50) rrhigh = 50; //2s intra beat interval -> 30 bpm
            if (rrmiss > 90) rrmiss = 90; //it means no beat for 3.6s (17bpm allowd, clearly out of range)
        }

        prevRegular = regular;
        if (rravg1 == rravg2)
        {
            regular = true;
        }
        // If the beat had been normal but turned odd, change the thresholds.
        else
        {
            regular = false;
            if (prevRegular)
            {
                threshold_i1 /= 2;
                threshold_f1 /= 2;
            }
        }
    }
   
    // At this point a back search can be implemented, in the case of a missing peak detection. (Drawback: The outputs must be delayed if this part is included)
    //For now it is not reported
    

    //Generating output signal
    if (qrs) outputSignal[current]=1;
    else outputSignal[current]=0;

    //Plotting output signal and averages
    int interval;
   
    // Giving as output rravg2 id the new RRinterval is higher than rrmiss. For now this feature is not implemented
    /*
    if(sample < 8){
        if (rr1[sample] != 0){
            interval = rr1[sample];
            HR = (600.0 * (float)(FS))/((float)(rr1[sample]));
        } else {
            HR = 0;
            interval =0;
        }
    
    } else{
        if (rr1[7] > rrmiss){
            if (rravg2 != 0){
                interval = rravg2;
                HR = (600.0 * (float)(FS))/((float)(rravg2));
            } else{
                HR = 0;
                interval = 0;
            }
        } else{
            if (rr1[7] != 0){
                interval = rr1[7];
                HR = (600.0 * (float)(FS))/((float)(rr1[7]));
            } else{
                interval =0;
                HR = 0;
            } 
        }
    }
    */
    if(sample < 8){
        if (rr1[sample] != 0){
            interval = rr1[sample];
            HR = (600.0 * (float)(FS))/((float)(rr1[sample]));
        } else {
            HR = 0;
            interval =0;
        }
    
    } else{
        
        if (rr1[7] != 0){
            interval = rr1[7];
            HR = (600.0 * (float)(FS))/((float)(rr1[7]));
        } else{
            interval =0;
            HR = 0;
        } 
        
    }

    output(outputSignal[current],interval,rr1[7],rravg1,rravg2);

    //For the median filter
    //Saving the results in a buffer to implement a moving average or median
    if(qrs){
        pt_results[pt_results_index] = HR; 
        pt_results_index++;
        if (pt_results_index >= HRM_HIST_LEN)
        {
            pt_results_index = 0;
        }
    }
    //Median
    //First: Sorting the buffer
    bool busy;
    do {
        busy = false;
        for (int i=0;i<HRM_HIST_LEN-1;i++) {
            if (pt_results[i] > pt_results[i+1]) {
                int te = pt_results[i];
                pt_results[i] = pt_results[i+1];
                pt_results[i+1] = te;
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
        if (pt_results[i]==0) continue;
        sumBPM += pt_results[i];
        n++;
    }
    if (n) {
        HR = (int)(sumBPM/n);
    } 
    
    return HR;



}

int output(ppg_t output, int interval, int interval_raw,int rravg1, int rravg2){
    #ifdef DUMP_FILE
        if (qrsFile)
        {
            if (!fprintf(qrsFile, "%d,%d,%d,%d,%d\n", output,interval,interval_raw,rravg1,rravg2))
                puts("error writing file");
            fflush(qrsFile);
        }
    #endif

}

