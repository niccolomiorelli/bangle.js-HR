#ifndef SPECTRALTRACKING_HEARTRATE_H
#define SPECTRALTRACKING_HEARTRATE_H

#include "../../types.h"


void spectralTracking_heartrate_init();

int spectralTracking_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int algorithm(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

#endif