#ifndef DEBUG1_HEARTRATE_H
#define DEBUG1_HEARTRATE_H

#include "../../types.h"

void debug1_heartrate_init();

int debug1_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_debug1(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);


#endif