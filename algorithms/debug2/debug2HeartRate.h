#ifndef DEBUG2_HEARTRATE_H
#define DEBUG2_HEARTRATE_H

#include "../../types.h"

void debug2_heartrate_init();

int debug2_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_debug2(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);


#endif