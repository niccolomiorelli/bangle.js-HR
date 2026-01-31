#ifndef DEBUG82_HEARTRATE_H
#define DEBUG82_HEARTRATE_H

#include "../../types.h"

void debug82_heartrate_init();

int debug82_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_debug82(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

#endif