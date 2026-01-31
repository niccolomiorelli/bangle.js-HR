#ifndef DEBUG10_HEARTRATE_H
#define DEBUG10_HEARTRATE_H

#include "../../types.h"

void debug10_heartrate_init();

int debug10_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_debug10(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

#endif