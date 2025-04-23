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
#include "preProcessingStage.h"
#include "config.h"

#include "../../utils/bandpass_filter/bandpass_filter.h"

//To print files
//#define DUMP_FILE
#ifdef DUMP_FILE
#include <stdio.h>
static FILE *magnitudeFile;
//static FILE *interpolatedFile;
#endif

static ring_buffer_t *inBuff;
static ring_buffer_t *outBuff;
static void (*nextStage)(void);
static uint8_t samplingPeriod = 40;    // in ms, this can be smaller than the actual sampling frequency, but it will result in more computations -> 40 ms
static uint16_t timeScalingFactor = 1; // use this for adjusting time to ms, in case the clock has higher precision
static time_delta_ms_t lastSampleTime = -1;

//Band Pass Filter
static BPFilter bpFilter;

// Variabile di stato per tracciare il tempo cumulativo
static uint64_t cumulative_time = 0;

void initPreProcessStage(ring_buffer_t *pInBuff, ring_buffer_t *pOutBuff, void (*pNextStage)(void))
{
    inBuff = pInBuff;
    outBuff = pOutBuff;
    nextStage = pNextStage;

    // Initialize the filter 
    BPFilter_init(&bpFilter);

    //Cumulative time to zero 
    cumulative_time = 0;

#ifdef DUMP_FILE
    magnitudeFile = fopen(DUMP_MAGNITUDE_FILE_NAME, "w+");
#endif
}

static void outPutDataPoint(data_point_t dp)
{
    lastSampleTime = dp.time;
    ring_buffer_queue(outBuff, dp);
    (*nextStage)();

}



void preProcessSample(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz)
{
    //FILTERING
    // Simply use the filter implemented in the bandpass.h and bandpass.c file
    BPFilter_put(&bpFilter, ppg);
    ppg_t ppg_filtered = BPFilter_get(&bpFilter);

    // Aggiornamento del tempo cumulativo
    // cumulative_time += (delta_ms / timeScalingFactor); //It was this before
    cumulative_time += (uint64_t)(delta_ms);

    // For now, this function doesnt do anything: Creating the dataPoint
    //What is lastSampleTime used for?

    data_point_t dataPoint;
    dataPoint.time = cumulative_time;
    dataPoint.magnitude = ppg_filtered;

#ifdef DUMP_FILE
    if (magnitudeFile)
    {
        if (!fprintf(magnitudeFile, "%lld, %d\n", dataPoint.time, dataPoint.magnitude))
            puts("error writing file");
    }
#endif

    outPutDataPoint(dataPoint);
}

void resetPreProcess(void)
{
    lastSampleTime = -1;

#ifdef DUMP_FILE
    if (magnitudeFile)
    {
        fflush(magnitudeFile);
    }
#endif
}