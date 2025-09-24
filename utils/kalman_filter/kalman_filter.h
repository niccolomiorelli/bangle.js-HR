#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#include "../../types.h"

double find_Qk(double acc_rms);
double find_Rk(double conf);

void kalman_HR_init(double* HR_est, double*P);

void kalman_HR_2meas_init(double* HR_est, double* P);

void kalman_HR_estimation(double HR_meas, double acc_rms, double conf, double* HR_est, double* P);

void kalman_HR_estimation_2measures(double HR_meas_2m_PPG, double HR_meas_2m_ACC, double acc_rms, double conf, double* HR_est_2m, double* P_2m);


#endif

