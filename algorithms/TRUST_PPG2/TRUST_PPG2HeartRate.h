#ifndef TRUST_PPG2_HEARTRATE_H
#define TRUST_PPG2_HEARTRATE_H

#include "../../types.h"

void trust_ppg2_heartrate_init();

int trust_ppg2_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

int main_algorithm_trust_ppg2(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);
#endif