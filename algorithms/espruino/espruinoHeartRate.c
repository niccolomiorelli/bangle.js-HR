/* ----------------------------------------------------
* ESPRUINO ALGORITHM
* ----------------------------------------------------
* This file is part of Espruino, a JavaScript interpreter for Microcontrollers
*
* Copyright (C) 2021 Gordon Williams <gw@pur3.co.uk>
*
* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, You can obtain one at http://mozilla.org/MPL/2.0/.
*
* ----------------------------------------------------------------------------
* Heart rate
* ----------------------------------------------------------------------------
*/

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <limits.h>

#include "espruinoHeartRate.h"
#include "../../types.h"


// Referring to the filter of order 175, not used for the moment
/*
 
 ==========================================================
 FIR filter designed with http://t-filter.engineerjs.com/
 
 AccelFilter_get modified to return 8 bits of fractional
 data.
 ==========================================================
 
 Source Code Tab:
 
 HRM
 Integer
 11
 int8_t
 int
 
 FIR filter designed with
  http://t-filter.appspot.com
 
 sampling frequency: 50 Hz
 
 fixed point precision: 11 bits
 
 * 0 Hz - 0.6 Hz
   gain = 0
   desired attenuation = -40 dB
   actual attenuation = n/a
 
 * 0.9 Hz - 3 Hz
   gain = 1
   desired ripple = 5 dB
   actual ripple = n/a
 
 * 3.6 Hz - 25 Hz
   gain = 0
   desired attenuation = -40 dB
   actual attenuation = n/a
 
*/

#define HRMFILTER_TAP_NUM 35 //175

 
typedef struct {
   ppg_t history[HRMFILTER_TAP_NUM];
   unsigned int last_index;
} HRMFilter;
 
/*
static const int8_t filter_taps[HRMFILTER_TAP_NUM] = {
   7,
   5,
   6,
   7,
   7,
   7,
   5,
   3,
   0,
   -3,
   -7,
   -10,
   -13,
   -15,
   -16,
   -15,
   -14,
   -11,
   -8,
   -4,
   0,
   3,
   5,
   6,
   6,
   5,
   3,
   1,
   -1,
   -2,
   -3,
   -2,
   -1,
   2,
   5,
   8,
   10,
   12,
   13,
   13,
   11,
   9,
   6,
   3,
   1,
   -1,
   -1,
   0,
   2,
   4,
   7,
   10,
   12,
   12,
   10,
   7,
   3,
   -2,
   -8,
   -13,
   -17,
   -19,
   -18,
   -16,
   -12,
   -8,
   -3,
   0,
   1,
   -1,
   -5,
   -13,
   -23,
   -33,
   -43,
   -51,
   -56,
   -55,
   -49,
   -37,
   -19,
   2,
   26,
   49,
   71,
   88,
   99,
   103,
   99,
   88,
   71,
   49,
   26,
   2,
   -19,
   -37,
   -49,
   -55,
   -56,
   -51,
   -43,
   -33,
   -23,
   -13,
   -5,
   -1,
   1,
   0,
   -3,
   -8,
   -12,
   -16,
   -18,
   -19,
   -17,
   -13,
   -8,
   -2,
   3,
   7,
   10,
   12,
   12,
   10,
   7,
   4,
   2,
   0,
   -1,
   -1,
   1,
   3,
   6,
   9,
   11,
   13,
   13,
   12,
   10,
   8,
   5,
   2,
   -1,
   -2,
   -3,
   -2,
   -1,
   1,
   3,
   5,
   6,
   6,
   5,
   3,
   0,
   -4,
   -8,
   -11,
   -14,
   -15,
   -16,
   -15,
   -13,
   -10,
   -7,
   -3,
   0,
   3,
   5,
   7,
   7,
   7,
   6,
   5,
   7
};
*/
static int16_t filter_taps[HRMFILTER_TAP_NUM]={
  -89, -72, -28,  2, -82, -336, -626, -661, -308, 104, -85, -1202,
-2609, -2858, -782, 3298, 7401, 9126, 7401, 3298, -782, -2858, -2609, -1202,
 -85, 104, -308, -661, -626, -336,  -82, 2, -28, -72, -89
};
 
static void HRMFilter_init(HRMFilter* f) {
   int i;
   for(i = 0; i < HRMFILTER_TAP_NUM; ++i)
     f->history[i] = 0;
   f->last_index = 0;
}
 
static void HRMFilter_put(HRMFilter* f, int input) {
   f->history[f->last_index++] = input;
   if(f->last_index == HRMFILTER_TAP_NUM)
     f->last_index = 0;
}
 
static int HRMFilter_get(HRMFilter* f) {
   long long acc = 0;
   int index = f->last_index, i;
   for(i = 0; i < HRMFILTER_TAP_NUM; ++i) {
     index = index != 0 ? index-1 : HRMFILTER_TAP_NUM-1;
     acc += (long long)f->history[index] * filter_taps[i];
   };
  //  int result = (int)(acc >> 4); //Nel caso del filtro precedente
   int result = (int)(acc >> 15);
   if (result > INT_MAX) result = INT_MAX; //I added this part
   if (result < INT_MIN) result = INT_MIN;
   return result;
}

HRMFilter hrmFilter;

///////////////////////////////////
#define HRM_HIST_LEN 16 
#define HRM_MEDIAN_LEN 8 

#define HRMVALUE_MIN -32768
#define HRMVALUE_MAX 32767


typedef struct {
  uint16_t bpm10; // 10x BPM
  uint8_t confidence; // 0..100%

  ppg_t raw;
  int16_t avg; // average signal value, moving average
  int16_t filtered;
  int16_t filtered1; // before filtered
  int16_t filtered2; // before filtered1
  bool wasLow; // has the signal gone below the average? set =false when a beat detected
  bool isBeat; // was this sample classified as a detected beat?
  int16_t times[HRM_HIST_LEN]; // times of previous beats, in 1/100th secs
  uint8_t timeIdx; // index in times
  int16_t timeBeat;
} HrmInfo;

int result = 0;

int compute_bpm(ppg_t x){
  return (int)(x*10);
}


HrmInfo hrmInfo = {0};


/// Initialise heart rate monitoring
void espruino_heartrate_init() {
  memset(&hrmInfo, 0, sizeof(hrmInfo));
  hrmInfo.wasLow = false;
  HRMFilter_init(&hrmFilter);
}

int16_t hrm_time_to_bpm10(int16_t time) {
  return (10 * 60 * 100) / time; // 10x BPM
}


bool hrm_had_beat() { 
  // Get time since last beat
  int16_t beatTime = (int16_t)(hrmInfo.timeBeat / 10); // in 1/100th sec
  hrmInfo.timeBeat = 0; //Riinizzializzo a zero
  if (beatTime<20) return false; // 1/5th sec is too short
  if (beatTime>255) beatTime=255;

  // store HRM times in list (in 1/100th sec)
  hrmInfo.times[hrmInfo.timeIdx] = (int16_t)beatTime;
  hrmInfo.timeIdx++;
  if (hrmInfo.timeIdx >= HRM_HIST_LEN)
    hrmInfo.timeIdx = 0;
  // copy times over
  int16_t times[HRM_HIST_LEN];
  memcpy(times, hrmInfo.times, sizeof(hrmInfo.times));
  // bubble sort
  bool busy;
  do {
    busy = false;
    for (int i=0;i<HRM_HIST_LEN-1;i++) {
      if (times[i] > times[i+1]) {
        int16_t t = times[i];
        times[i] = times[i+1];
        times[i+1] = t;
        busy = true;
      }
    }
  } while (busy);
  // calculate HRM from middle values
  int16_t min = (int16_t)((HRM_HIST_LEN - HRM_MEDIAN_LEN)/2); // 4
  int16_t max = (int16_t)((HRM_HIST_LEN + HRM_MEDIAN_LEN)/2); // 12
  int16_t n = 0;
  int16_t sumBPM = 0;
  for (int16_t i=min;i<max;i++) {
    if (times[i]==0) continue;
    int16_t BPM10 = hrm_time_to_bpm10(times[i]); // 10x BPM
    sumBPM += BPM10;
    n++;
  }
  if (n) {
    hrmInfo.bpm10 = (int16_t)(sumBPM/n);
    if (n >= HRM_MEDIAN_LEN) { // not enough values to be confident
      // spread = difference between min+max BPM*10
      int spread = hrm_time_to_bpm10(times[min]) - hrm_time_to_bpm10(times[max]);
      if (spread > 100) spread -= 100; // 10bpm difference = 100% accuracy
      else spread = 0;
      if (spread > 400) spread=400; // 40bpm difference = low accuracy
      hrmInfo.confidence = 100 - spread/4;
      if (hrmInfo.bpm10 < 300)
        hrmInfo.confidence = 0; // not confident about BPM less than 30!
    } else
      hrmInfo.confidence = 0;
  } else {
    hrmInfo.confidence = 0;
  }
  if (n!=0) { // do we have a useful HRM value?
   return true; //It means the bpm10 value has been updated
   } else return false;

}

/// Add new heart rate value
int espruino_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz) {
  if (ppg<HRMVALUE_MIN) ppg=HRMVALUE_MIN;
  if (ppg>HRMVALUE_MAX) ppg=HRMVALUE_MAX;
  hrmInfo.raw = ppg;
  HRMFilter_put(&hrmFilter, ppg);
  int h = HRMFilter_get(&hrmFilter);
  if (h<=-32768) h=-32768;
  if (h>32767) h=32767;
  hrmInfo.filtered2 = hrmInfo.filtered1;
  hrmInfo.filtered1 = hrmInfo.filtered;  
  hrmInfo.filtered = h;

  //
  bool hadBeat = false;
  hrmInfo.isBeat = false;

  // Cumulative time
  hrmInfo.timeBeat += delta_ms;

  if (h < hrmInfo.avg)
    hrmInfo.wasLow = true;
  else if (hrmInfo.wasLow && (hrmInfo.filtered1 >= hrmInfo.filtered) && (hrmInfo.filtered1 >= hrmInfo.filtered2)) {
    hrmInfo.wasLow = false; // peak detected, and had previously gone below average
    hrmInfo.isBeat = true;
    hadBeat = hrm_had_beat();
  }
  


  hrmInfo.avg = ((hrmInfo.avg*7) + h) >> 3;

  return hrmInfo.bpm10;
  
}






