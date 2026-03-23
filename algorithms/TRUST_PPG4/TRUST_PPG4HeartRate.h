#ifndef TRUST_PPG4_HEARTRATE_H
#define TRUST_PPG4_HEARTRATE_H

#include "../../types.h"

void trust_ppg4_heartrate_init();

int trust_ppg4_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_trust_ppg4(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

#endif