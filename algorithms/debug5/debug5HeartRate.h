#ifndef DEBUG5_HEARTRATE_H
#define DEBUG5_HEARTRATE_H

#include "../../types.h"

void debug5_heartrate_init();

int debug5_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_debug5(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);


#endif