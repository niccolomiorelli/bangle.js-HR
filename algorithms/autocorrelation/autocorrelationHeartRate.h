#ifndef AUTOCORRELATION_HEARTRATE_H
#define AUTOCORRELATION_HEARTRATE_H

#include "../../types.h"



 /// Initialise heart rate monitoring
 void autocorrelation_heartrate_init();

 /// Add new heart rate value, return true if there was a heart beat
 int autocorrelation_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);
 
 #endif