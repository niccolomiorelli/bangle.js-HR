#ifndef DEBUG6_HEARTRATE_H
#define DEBUG6_HEARTRATE_H

#include "../../types.h"

void debug6_heartrate_init();

int debug6_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_debug6(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);


#endif