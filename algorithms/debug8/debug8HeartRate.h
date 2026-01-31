#ifndef DEBUG8_HEARTRATE_H
#define DEBUG8_HEARTRATE_H

#include "../../types.h"

void debug8_heartrate_init();

int debug8_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_debug8(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

#endif