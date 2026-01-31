#ifndef DEBUG892_HEARTRATE_H
#define DEBUG892_HEARTRATE_H

#include "../../types.h"

void debug892_heartrate_init();

int debug892_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_debug892(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);
#endif