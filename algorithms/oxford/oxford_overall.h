/*
The MIT License (MIT)

Copyright (c) 2020 Anna Brondin and Marcus Nordström

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef STEP_COUNTING_ALGO_H
#define STEP_COUNTING_ALGO_H
#include <stdint.h>
#include "../../types.h"
#include "config.h"

/**
    Initializes all buffers and everything the algorithm needs
*/
void oxford_init(void);

/**
    Resets the Heart Rate (For now it is taken from the step count which makes sense, I dont know yet if it makes sense with HR)
*/
void oxford_resetHR(void);

/**
    Resets the entire algorithm
*/
void oxford_resetAlgo(void);

/**
    This function takes the raw accelerometry data and computes the entire algorithm
    @param time, the current time in ms
    @param x, the x axis
    @param y, the y axis
    @param z, the z axis
    @return steps walked
*/
int oxford_HR(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);

#endif
