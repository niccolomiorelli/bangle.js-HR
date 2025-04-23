/*
The MIT License (MIT)

Copyright (c) 2020 Anna Brondin and Marcus Nordström

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>

#include "oxford_overall.h"
#include "ringbuffer.h"
#include "preProcessingStage.h"
#include "scoringStage.h"
#include "detectionStage.h"
#include "postProcessingStage.h"

// General data
static int HR; //THIS WOULD BE oxford_HR (the variable, not the function)
// Buffers
static ring_buffer_t rawBuf;
static ring_buffer_t ppBuf;
static ring_buffer_t peakScoreBuf;
static ring_buffer_t peakBuf;

//For the median filter
//Buffer per implementare media mobile o mediana sui risultati
#define HRM_HIST_LEN 8
#define HRM_MEDIAN_LEN 4
static int fft_results[HRM_HIST_LEN] = {0};
static int fft_results_index = 0;

static void increaseStepCallback(uint64_t interval)
{
    // print for debugging if I need to
    HR = (int)(600000.0/(float)(interval)); //HR = 1/interval * 60 * 1000 (to convert from ms) * 10 to obtain bpm*10

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

void oxford_init()
{
    oxford_resetHR();
    // init buffers
    ring_buffer_init(&rawBuf);
    ring_buffer_init(&ppBuf);
    ring_buffer_init(&peakScoreBuf);
    ring_buffer_init(&peakBuf);

    initPreProcessStage(&rawBuf, &ppBuf, scoringStage);
    initScoringStage(&ppBuf, &peakScoreBuf, detectionStage);
    initDetectionStage(&peakScoreBuf, &peakBuf, postProcessingStage);
    initPostProcessingStage(&peakBuf, &increaseStepCallback);
    
}

void oxford_resetHR(void)
{
    HR = 0;
}

void oxford_resetAlgo(void)
{
    resetPreProcess();
    resetDetection();
    resetPostProcess();
    ring_buffer_init(&rawBuf);
    ring_buffer_init(&ppBuf);
    ring_buffer_init(&peakScoreBuf);
    ring_buffer_init(&peakBuf);
}

int oxford_HR(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz)
{
    preProcessSample(delta_ms, ppg, accx, accy, accz);
    return HR;
}