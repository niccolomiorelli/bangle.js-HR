#ifndef DEBUG11_HEARTRATE_H
#define DEBUG11_HEARTRATE_H

#include "../../types.h"

void debug11_heartrate_init();

int debug11_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_debug11(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);
#endif