#ifndef TRUST_PPG3_HEARTRATE_H
#define TRUST_PPG3_HEARTRATE_H

#include "../../types.h"

void trust_ppg3_heartrate_init();

int trust_ppg3_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_trust_ppg3(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

#endif