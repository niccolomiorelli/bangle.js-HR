/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * ----------------------------------------------------------------------------
 * Dummy Step Counter
 * ----------------------------------------------------------------------------
 */

#ifndef DUMMY_HEARTRATE_H
#define DUMMY_HEARTRATE_H

#include "../../types.h"

/**
 * Initalise and reset the step counter
 */
void dummy_heartrate_init();

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
int dummy_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

#endif
