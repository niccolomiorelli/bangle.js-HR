////////////////////////////////////////////////////////////////////////////////////////
/*
 * This file is part of Espruino, a JavaScript interpreter for Microcontrollers
 *
 * Copyright (C) 2021 Gordon Williams <gw@pur3.co.uk>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * ----------------------------------------------------------------------------
 * heart rate monitoring
 * ----------------------------------------------------------------------------
*/

#ifndef ESPRUINO_HEARTRATE_H
#define ESPRUINO_HEARTRATE_H

#include "../../types.h"
#include <time.h>
#include <stdbool.h>



 /// Initialise heart rate monitoring
 void espruino_heartrate_init();
 
 bool hrm_had_beat();

 /// Add new heart rate value, return true if there was a heart beat
 int espruino_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);
 
 #endif





