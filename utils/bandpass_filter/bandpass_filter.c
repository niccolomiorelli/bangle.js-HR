/* ----------------------------------------------------
* Bandpass Filter
* ----------------------------------------------------
* 
* Filter designed with the python ibrary "scipy.signal", the coefficients are extracted using the function firwin()
*
* Filter characteristics:
* Sampling freuency: 25 Hz
* Cut frequencies: (0.5, 4)Hz
* Order: 35
* Window type: Hamming
*
* Fixed point precision: 15 bits
* The coefficients are multiplied by 2^15, so then the filtered value (resulting dìfrom the convolution of the signal with the coefficients) are divided by 2^15
*
*/

#include "bandpass_filter.h"
#include <limits.h>

static int16_t filter_taps[BPFILTER_TAP_NUM]={
    -89, -72, -28,  2, -82, -336, -626, -661, -308, 104, -85, -1202,
 -2609, -2858, -782, 3298, 7401, 9126, 7401, 3298, -782, -2858, -2609, -1202,
   -85, 104, -308, -661, -626, -336,  -82, 2, -28, -72, -89
};

//Filter not used for now:
/*
static int16_t filter_taps[BPFILTER_TAP_NUM] = {
    -13,
    3,
    -73,
    -214,
    -254,   
    -94,   
    -10,  
    -410, 
    -1068,
    -1086, 
    -181,   
    345,
    -1000, 
    -3282, 
    -3089,  
    1710,  
    8530, 
    11775,  
    8530,  
    1710, 
    -3089, 
    -3282, 
    -1000,   
    345,
    -181, 
    -1086, 
    -1068,  
    -410,   
    -10,   
    -94,  
    -254,  
    -214,   
    -73,     
    3,   
    -13};
*/

void BPFilter_init(BPFilter *f)
{
    int i;
    for (i = 0; i < BPFILTER_TAP_NUM; ++i)
        f->history[i] = 0;
    f->last_index = 0;
}

void BPFilter_put(BPFilter *f, ppg_t input)
{
    f->history[f->last_index++] = input;
    if (f->last_index == BPFILTER_TAP_NUM)
        f->last_index = 0;
}

int BPFilter_get(BPFilter *f)
{
    long long acc = 0;
    int index = f->last_index, i;
    for (i = 0; i < BPFILTER_TAP_NUM; ++i)
    {
        index = index != 0 ? index - 1 : BPFILTER_TAP_NUM - 1;
        acc += (long long)f->history[index] * filter_taps[i];
    };
    int temp = acc >> 15;  //NB: Shift of 15 because to create the filter coefficients I have moltiplied by 2^15
    if (temp > INT_MAX) temp = INT_MAX; //I added this part to avoid overflow
    if (temp < INT_MIN) temp = INT_MIN;
    return temp;
}