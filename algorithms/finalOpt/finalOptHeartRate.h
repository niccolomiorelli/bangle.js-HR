#ifndef FINALOPT_HEARTRATE_H
#define FINALOPT_HEARTRATE_H

#include "../../types.h"

void finalOpt_heartrate_init();

int finalOpt_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_opt(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);


#endif