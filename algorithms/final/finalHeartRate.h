#ifndef FINAL_HEARTRATE_H
#define FINAL_HEARTRATE_H

#include "../../types.h"

void final_heartrate_init();

int final_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);


#endif