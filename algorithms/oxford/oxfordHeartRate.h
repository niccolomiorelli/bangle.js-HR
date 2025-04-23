/**
 * Wrapper of the Oxford step counter.
 */
#ifndef OXFORD
#define OXFORD

#include "./oxford_overall.h"

// Initialise step counting
void oxford_heartrate_init()
{
    oxford_init();
    oxford_resetHR();
    oxford_resetAlgo();
}

// process sample
int oxford_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz)
{
    return oxford_HR(delta_ms, ppg, accx, accy, accz);
}

#endif