/* ----------------------------------------------------
* Kalman Filter
* ----------------------------------------------------
* 

*/

#include "kalman_filter.h"
#include <math.h>

#define DELTA_MIN 0.05    //Hz
#define DELTA_MAX 0.83    //Hz
#define Z_SCORE   1.96    //97.5% confidence interval
#define MIN_CONF  0.001   //Safety threshold for confidence
#define THRESHOLD1_QK 8500.0
#define THRESHOLD2_QK 11000.0

//Variance of the state: how much HR oscillates: it is a function of the acc (rms_acc) -> STEP-WISE function in this case
double find_Qk(double acc_rms)
{
    double var = 0.0;

    if (acc_rms <= THRESHOLD1_QK) //Riposo o camminata leggera
        var = 0.02; //in Hz, so 1.2 bpm (referring to 2.5s windows)
    else if (acc_rms <= THRESHOLD2_QK) //Camminata moderata
        var = 0.05; //#in Hz, so 3bom
    else
        var = 0.08; //in Hz, so 4.8 bpm

    return var * var; //Qk = variance^2
}

//The idea is: HR_est +- delta = HR_est +- 1.96*sigma, assuming a normal distribution: 1.96 represents the 97.5% percentile 
double find_Rk(double conf)
{
    if (conf < MIN_CONF)
        conf = MIN_CONF;

    double delta_Hz = DELTA_MIN + (1.0f - conf) * (DELTA_MAX - DELTA_MIN);
    double sigma = delta_Hz / Z_SCORE;
    return (sigma * sigma); //Rk = sigma^2
}

void kalman_HR_init(double* HR_est, double* P)
{
    *HR_est = 1.0;
    *P = 0.2;
}

void kalman_HR_2meas_init(double* HR_est, double* P)
{
    *HR_est = 1.0;
    *P = 0.1;
}

void kalman_HR_estimation(double HR_meas, double acc_rms, double conf, double* HR_est, double* P)
{
    double Qk = find_Qk(acc_rms);
    double Rk = find_Rk(conf);

    //Prediction
    double HR_pred = *HR_est;
    double P_pred = *P + Qk;

    if (conf > 0.3){
    //Update
    double K = P_pred / (P_pred + Rk);
    *HR_est = HR_pred + K * (HR_meas - HR_pred);
    *P = (1.0 - K) * P_pred;
    } else {
        *HR_est = HR_pred;
        *P = P_pred;
    }
}

void kalman_HR_estimation_2measures(double HR_meas_2m_PPG, double HR_meas_2m_ACC, double acc_rms, double conf, double* HR_est_2m, double* P_2m)
{
    
    double Qk  = find_Qk(acc_rms);
    double Rk_PPG = find_Rk(conf);
    double Rk_ACC = (0.25/1.96)*(0.25/1.96);

    //Prediction:
    double HR_pred = *HR_est_2m;
    double P_pred = *P_2m + Qk;

    if (conf >= 0.3){

    //Update
    double det  = P_pred*(Rk_PPG+Rk_ACC) + Rk_PPG*Rk_ACC;
    double K1 = P_pred*Rk_ACC / det;
    double K2 = P_pred*Rk_PPG / det;

    //Innovation
    double y1 = HR_meas_2m_PPG - HR_pred;
    double y2 = HR_meas_2m_ACC - HR_pred;

    *HR_est_2m = HR_pred + (K1 * y1 + K2 * y2);

    double KH = K1 + K2;
    *P_2m = (1.0 - KH) * P_pred;

    } else {
        *HR_est_2m = HR_pred;
        *P_2m = P_pred;
    }

}
