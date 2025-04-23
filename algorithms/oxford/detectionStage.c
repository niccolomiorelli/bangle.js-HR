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

#include "detectionStage.h"
#include "postProcessingStage.h"
#include "config.h"

//IMPORTANT: THRESHOLDS TO SET --change here--
#define THRESHOLD_DER 20 //For the condition on the steep derivative
#define THRESHOLD_STD 60 //For the condition on the std, in order to use the fast or slow alpha



#ifdef DUMP_FILE
#include <stdio.h>
static FILE
    *detectionFile;
static FILE *detectionFile;
static FILE 
    *meanFile;
static FILE *meanFile;
#endif

static ring_buffer_t *inBuff;
static ring_buffer_t *outBuff;
static void (*nextStage)(void);

//Adding a new buffer for the second condition: Steep derivative 
static ring_buffer_t derBuffer;

static float mean = 0.0;
static float std = 0.0;
static float variance = 0.0;
static int16_t count = 0;
static int16_t threshold_int = DETECTION_TRHE_WHOLE;
static int16_t threshold_frac = DETECTION_TRHE_PART;

//For EMA (initialize with the slow alpha)
float alpha = 0.03;

void initDetectionStage(ring_buffer_t *pInBuff, ring_buffer_t *peakBufIn, void (*pNextStage)(void))
{
    inBuff = pInBuff;
    outBuff = peakBufIn;
    nextStage = pNextStage;
    count = 0;
    alpha = 0.03; //Initialize to the slow alpha
    //For the derivative condition:
    ring_buffer_init(&derBuffer);
    
#ifdef DUMP_FILE
    detectionFile = fopen(DUMP_DETECTION_FILE_NAME, "w+");
    meanFile = fopen(DUMP_MEAN_FILE_NAME, "w+");
#endif
}

void detectionStage(void)
{
    if (!ring_buffer_is_empty(inBuff))
    {
        //NOTE: the mean and std are used as float numbers, to avoid problems during their computation
        float oMean = mean;
        data_point_t dataPoint;
        ring_buffer_dequeue(inBuff, &dataPoint);
        count++;

        //Updating the derBuffer
        ring_buffer_queue(&derBuffer,dataPoint);

        // APPROACH 1: Welford method: Compute the mean and std of the entire signal
        /*
        if (count == 1)
        {
            mean = (float)dataPoint.magnitude;
            std = 0;
        }
        else if (count == 2)
        {
            mean = (mean + (float)dataPoint.magnitude) / 2.0;
            std = (float)sqrt((((float)dataPoint.magnitude - mean) * ((float)dataPoint.magnitude - mean)) + ((oMean - mean) * (oMean - mean))) / 2.0;
        }
        else
        {
            mean = ((float)dataPoint.magnitude + (((float)count - 1) * mean)) / ((float)count);
            //mean = dataPoint.magnitude;
            // ppg_t part1 = ((std * std) / (count - 1)) * (count - 2);
            // ppg_t part2 = ((oMean - mean) * (oMean - mean));
            // ppg_t part3 = ((dataPoint.magnitude - mean) * (dataPoint.magnitude - mean)) / count;
            float part1 = ((std * std) / (count - 1)) * (count - 2);
            float part2 = ((oMean - mean) * (oMean - mean));
            float part3 = ((dataPoint.magnitude - mean) * (dataPoint.magnitude - mean)) / count;
            std = (float)sqrt(part1 + part2 + part3);
        }
        */

        // APPROACH 2: EMA: Exponential Moving Average
        //Using two alpha values: one slow and one fast. In this way, noisy portions of the signal impact less on the mean and std and the ricovery is faster
        if (std > THRESHOLD_STD){ 
            alpha = 0.3;
        } else {
            alpha = 0.03;
        }

        if (count ==1){
            mean = (float)dataPoint.magnitude;
            std = 0;
        }
        else{
            mean = alpha*(float)dataPoint.magnitude + (1-alpha)*mean;
            variance = alpha*((float)dataPoint.magnitude - mean)*((float)dataPoint.magnitude -mean) + (1-alpha)*variance;
            std = sqrt(variance);
        }
        

#ifdef DUMP_FILE
        if (meanFile)
        {
            // Printing the mean and std also
            if (!fprintf(meanFile, "%lld, %d, %f, %f\n", dataPoint.time, dataPoint.magnitude,mean,std))
                puts("error writing file");
            fflush(detectionFile);
        }
#endif
        //To use the positive threshold
        float threshold_float = std*0.9;
        //To use the negative threshold with THE REVERESE LOGIC
        //float threshold_float = -std*1.3;
        if (count > 15)
        {
            if (((float)dataPoint.magnitude - mean) < threshold_float)
            {
                //I have to do another conditioin: for the steep derivative
                data_point_t dataPoint2;
                uint8_t numDer = ring_buffer_num_items(&derBuffer);
                //Fetching the three last elements in derBuffer
                ppg_t lasts[4];
                for (uint8_t i=0;i<4;i++){
                    ring_buffer_peek(&derBuffer,&dataPoint2,numDer-1-i);
                    lasts[i] = dataPoint2.magnitude; //where lasts[0] is the last element, lasts[1] is the second last and lasts[3] is the third last and so on
                }

                ppg_t der=0;
                if(lasts[3]<=lasts[2] && lasts[2]<=lasts[1] && lasts[1]<=lasts[0]){
                    der = lasts[0] - lasts[3];
                }else if (lasts[3]>=lasts[2] && lasts[2]<=lasts[1] && lasts[1]<=lasts[0]){
                    der = lasts[0] - lasts[2];
                }else if (lasts[2]>=lasts[1] && lasts[1]<=lasts[0]){   
                    der = lasts[0] - lasts[1];
                }
            
                //POSITIVE THRESHOLD
                //Second condition on the derivative:
                if(der > THRESHOLD_DER){
                    // This is a peak:
                    ring_buffer_queue(outBuff, dataPoint);
                    (*nextStage)();
                }
                

               //NRGATIVE THRESHOLD
               /*
               ring_buffer_queue(outBuff, dataPoint);
               (*nextStage)();
               */

#ifdef DUMP_FILE
                if (detectionFile)
                {
                    if (!fprintf(detectionFile, "%lld, %d, %f, %f,%f,%d,%d,%d,%d,%d,%d\n", dataPoint.time, dataPoint.magnitude, threshold_float,mean,std,numDer,lasts[3],lasts[2],lasts[1],lasts[0],der))
                        puts("error writing file");
                    fflush(detectionFile);
                }
#endif
            }
        }
    }
}

void resetDetection(void)
{
    std = 0;
    mean = 0;
    variance = 0;
    count = 0;
}
