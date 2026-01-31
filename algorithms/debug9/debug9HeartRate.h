#ifndef DEBUG9_HEARTRATE_H
#define DEBUG9_HEARTRATE_H

#include "../../types.h"

void debug9_heartrate_init();

int debug9_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_debug9(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

#endif