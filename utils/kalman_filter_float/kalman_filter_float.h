#ifndef KALMAN_FILTER_FLOAT_H_
#define KALMAN_FILTER_FLOAT_H_

#include "../../types.h"

float find_Qk_float(float acc_rms);
float find_Rk_float(float conf);

void kalman_HR_init_float(float* HR_est, float*P);

void kalman_HR_2meas_init_float(float* HR_est, float* P);

void kalman_HR_estimation_float(float HR_meas, float acc_rms, float conf, float* HR_est, float* P);

void kalman_HR_estimation_2measures_float(float HR_meas_2m_PPG, float HR_meas_2m_ACC, float acc_rms, float conf, float* HR_est_2m, float* P_2m);


#endif

