#ifndef DEBUG88_HEARTRATE_H
#define DEBUG88_HEARTRATE_H

#include "../../types.h"

void debug88_heartrate_init();

int debug88_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_debug88(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

#endif