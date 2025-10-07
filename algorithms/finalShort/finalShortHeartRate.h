#ifndef FINALSHORT_HEARTRATE_H
#define FINALSHORT_HEARTRATE_H

#include "../../types.h"

void finalShort_heartrate_init();

int finalShort_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

#endif
