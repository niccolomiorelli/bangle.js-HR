#ifndef PANTOMPKINS_HEARTRATE_H
#define PANTOMPKINS_HEARTRATE_H

#include "../../types.h"

/**
 * Initalise and reset the step counter
 */
void pantompkins_heartrate_init();

/* Registers a new data point for heart rate.
 *
 * delta_ms: difference in millis between this sample and the previous
 * ppg: photoplethysmogram value
 * accx: acceleration on the x axis
 * accy: acceleration on the y axis
 * accz: acceleration on the z axis
 *
 * returns: the number of steps counted so far.
 */
int pantompkins_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

#endif