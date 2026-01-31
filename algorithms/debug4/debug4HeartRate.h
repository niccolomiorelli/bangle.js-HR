#ifndef DEBUG4_HEARTRATE_H
#define DEBUG4_HEARTRATE_H

#include "../../types.h"

void debug4_heartrate_init();

int debug4_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_debug4(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);


#endif