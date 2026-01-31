#ifndef DEBUG89_HEARTRATE_H
#define DEBUG89_HEARTRATE_H

#include "../../types.h"

void debug89_heartrate_init();

int debug89_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_debug89(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);
#endif