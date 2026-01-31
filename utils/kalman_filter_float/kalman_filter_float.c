/* ----------------------------------------------------
* Kalman Filter FLOAT version
* ----------------------------------------------------
* 

*/

#include "kalman_filter_float.h"
#include <math.h>

#define DELTA_MIN 0.05f    //Hz
#define DELTA_MAX 0.83f    //Hz
#define Z_SCORE   1.96f    //97.5% confidence interval
#define MIN_CONF  0.001f   //Safety threshold for confidence
#define THRESHOLD1_QK 8500.0f
#define THRESHOLD2_QK 11000.0f
#define MAX_DHR_HZ 0.05f
#define STATE_STATIONARY 0
#define STATE_TRANSITION 1
#define STATE_MOTION_LOW 2
#define STATE_MOTION_MEDIUM 3
#define STATE_MOTION_HIGH 4

float clampf(float x, float lo, float hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}

//Variance of the state: how much HR oscillates: it is a function of the acc (rms_acc) -> STEP-WISE function in this case
float find_Qk_float(float acc_rms)
{
    float var = 0.0f;

    if (acc_rms <= THRESHOLD1_QK) //Riposo o camminata leggera
        var = 0.02f; //in Hz, so 1.2 bpm (referring to 2.5s windows)
    else if (acc_rms <= THRESHOLD2_QK) //Camminata moderata
        var = 0.05f; //#in Hz, so 3bom
    else
        var = 0.08f; //in Hz, so 4.8 bpm

    return var * var; //Qk = variance^2
}
float find_Qk_float2(int state)
{
    float var = 0.0f;

    if (state  == STATE_STATIONARY || state == STATE_TRANSITION || state == STATE_MOTION_LOW) //Riposo o camminata leggera
        var = 0.01f; //in Hz, so  bpm (referring to 2.5s windows)
    else if (state == STATE_MOTION_MEDIUM) //Camminata moderata
        var = 0.02f; //#in Hz, so bom
    else
        var = 0.03f; //in Hz, so bpm

    return var * var; //Qk = variance^2
}

//The idea is: HR_est +- delta = HR_est +- 1.96*sigma, assuming a normal distribution: 1.96 represents the 97.5% percentile 
float find_Rk_float(float conf)
{
    if (conf < MIN_CONF)
        conf = MIN_CONF;

    float delta_Hz = DELTA_MIN + (1.0f - conf) * (DELTA_MAX - DELTA_MIN);
    float sigma = delta_Hz / Z_SCORE;
    return (sigma * sigma); //Rk = sigma^2
}
float find_Rk_float2(float conf)
{
    if (conf < MIN_CONF)
        conf = MIN_CONF;

    float delta_Hz = DELTA_MIN + (1.0f - conf) * (1.2f - DELTA_MIN);
    float sigma = delta_Hz / Z_SCORE;
    return (sigma * sigma); //Rk = sigma^2
}

void kalman_HR_init_float(float* HR_est, float* P)
{
    *HR_est = 1.0f;
    *P = 0.2f;
}

void kalman_HR_2meas_init_float(float* HR_est, float* P)
{
    *HR_est = 1.0f;
    *P = 0.1f;
}

void kalman_HR_estimation_float(float HR_meas, float acc_rms, float conf, float* HR_est, float* P)
{
    float Qk = find_Qk_float(acc_rms);
    float Rk = find_Rk_float(conf);

    //Prediction
    float HR_pred = *HR_est;
    float P_pred = *P + Qk;

    if (conf > 0.3f){
    //Update
    float K = P_pred / (P_pred + Rk);
    *HR_est = HR_pred + K * (HR_meas - HR_pred);
    *P = (1.0f - K) * P_pred;
    } else {
        *HR_est = HR_pred;
        *P = P_pred;
    }
}
void kalman_HR_estimation_float2(float HR_meas, int state, float conf, float* HR_est, float* P)
{
    float Qk = find_Qk_float2(state);
    float Rk = find_Rk_float2(conf);

    //Prediction
    float HR_old = *HR_est;


    float HR_pred = *HR_est;
    float P_pred = *P + Qk;

    if (conf > 0.3f){
    //Update
    float K = P_pred / (P_pred + Rk);
    *HR_est = HR_pred + K * (HR_meas - HR_pred);
    *P = (1.0f - K) * P_pred;
    } else {
        *HR_est = HR_pred;
        *P = P_pred;
    }

    // // HARD boundary: limit step-to-step change
    // float d = *HR_est - HR_old;
    // d = clampf(d, -MAX_DHR_HZ, +MAX_DHR_HZ);
    // *HR_est = HR_old + d;
}

void kalman_HR_estimation_2measures_float(float HR_meas_2m_PPG, float HR_meas_2m_ACC, float rms_acc, float conf, float* HR_est_2m, float* P_2m)
{
    
    float Qk  = find_Qk_float(rms_acc);
    float Rk_PPG = find_Rk_float(conf);
    float Rk_ACC = (0.25f/1.96f)*(0.25f/1.96f);

    //Prediction:
    float HR_pred = *HR_est_2m;
    float P_pred = *P_2m + Qk;

    if (conf >= 0.3){

        //Update
        float det  = P_pred*(Rk_PPG+Rk_ACC) + Rk_PPG*Rk_ACC;
        float K1 = P_pred*Rk_ACC / det;
        float K2 = P_pred*Rk_PPG / det;

        //Innovation
        float y1 = HR_meas_2m_PPG - HR_pred;
        float y2 = HR_meas_2m_ACC - HR_pred;

        *HR_est_2m = HR_pred + (K1 * y1 + K2 * y2);

        float KH = K1 + K2;
        *P_2m = (1.0f - KH) * P_pred;

    } else {
        *HR_est_2m = HR_pred;
        *P_2m = P_pred;
    }

}
void kalman_HR_estimation_2measures_float2(float HR_meas_2m_PPG, float HR_meas_2m_ACC, int state, float conf, float* HR_est_2m, float* P_2m)
{
    
    float Qk  = find_Qk_float2(state);
    float Rk_PPG = find_Rk_float(conf);
    float Rk_ACC = (0.25f/1.96f)*(0.25f/1.96f);

    float HR_old_2m = *HR_est_2m;

    //Prediction:
    float HR_pred = *HR_est_2m;
    float P_pred = *P_2m + Qk;

    if (conf >= 0.3){

        //Update
        float det  = P_pred*(Rk_PPG+Rk_ACC) + Rk_PPG*Rk_ACC;
        float K1 = P_pred*Rk_ACC / det;
        float K2 = P_pred*Rk_PPG / det;

        //Innovation
        float y1 = HR_meas_2m_PPG - HR_pred;
        float y2 = HR_meas_2m_ACC - HR_pred;

        *HR_est_2m = HR_pred + (K1 * y1 + K2 * y2);

        float KH = K1 + K2;
        *P_2m = (1.0f - KH) * P_pred;

    } else {
        *HR_est_2m = HR_pred;
        *P_2m = P_pred;
    }

    // HARD boundary: limit step-to-step change
    // float d = *HR_est_2m - HR_old_2m;
    // d = clampf(d, -MAX_DHR_HZ, +MAX_DHR_HZ);
    // *HR_est_2m = HR_old_2m  + d;

}