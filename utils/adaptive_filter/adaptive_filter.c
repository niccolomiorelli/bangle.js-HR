/* ----------------------------------------------------
* Adaptive Filter
* ----------------------------------------------------
* 

*/

#include "adaptive_filter.h"
#include <limits.h>

void AdaptFilter_init(AdaptFilter *f)
{
    int i;
    for (i = 0; i < ORDER; ++i){
        f->x[i] = 0.0;
        f->weights[i] = 0.0;
    }
    f->d = 0.0;
    f->last_index = 0;
}

void AdaptFilter_put(AdaptFilter *f, float x, float d)
{
    int i = f->last_index;
    f->x[i] = x;
    f->last_index = (i + 1) % ORDER;
    f->d = d;
}

/*
float AdaptFilter_get(AdaptFilter *f)
{
    float xvec[ORDER];
    int index = f->last_index;
    for (int i = 0; i < ORDER; ++i) {
        index = (index != 0) ? index - 1 : ORDER - 1;
        xvec[i] = f->x[index]; //Buffer with order samples
    }

    float acc = 0;
    for (int i = 0; i < ORDER; ++i) {
        acc += f->weights[i] * xvec[i];
    }

    //Compute the output
    float e = f->d - acc;
    float mu = 0.001;

    // Aggiorna i pesi
    for (int i = 0; i < ORDER; ++i) {
        f->weights[i] += mu * e * xvec[i];
    }

    return e;
}
*/
float AdaptFilter_get(AdaptFilter *f)
{
    float xvec[ORDER];
    int index = f->last_index;

    for (int i = 0; i < ORDER; ++i) {
        index = (index != 0) ? index - 1 : ORDER - 1;
        xvec[i] = f->x[index];
    }

    float acc = 0.0f;
    for (int i = 0; i < ORDER; ++i) {
        acc += f->weights[i] * xvec[i];
    }

    float e = f->d - acc;

    float mu = 0.0001f;         
    float epsilon = 1e-6f;   

    float norm = epsilon;
    for (int i = 0; i < ORDER; ++i) {
        norm += xvec[i] * xvec[i];
    }

    // Aggiorna i pesi usando NLMS
    for (int i = 0; i < ORDER; ++i) {
        f->weights[i] += (mu) * e * xvec[i];
    }

    return e;
}


/*
int AdaptFilter_get(AdaptFilter *f)
{
    long long acc = 0;
    int index = f->last_index;
    for (int i = 0; i < ORDER; ++i)
    {
        index = index != 0 ? index - 1 : ORDER - 1; //It means 
        acc += (long long)f->x[index] * f->weights[i];
    };
    int y_hat = acc >> 15;  //NB: Shift of 15 because to create the filter coefficients I have moltiplied by 2^15
    if (y_hat > INT_MAX) y_hat = INT_MAX; //I added this part to avoid overflow
    if (y_hat < INT_MIN) y_hat = INT_MIN;
    
    int16_t e = f->d - (int16_t)y_hat;


    for (int i = 0; i < ORDER; ++i) {
        index = (index != 0) ? index - 1 : ORDER - 1;
        // Ridimensionamento e aggiornamento del coefficiente
        f->weights[i] += (MU * e * f->x[index]) >> 15;
    }

    return e;

}
*/