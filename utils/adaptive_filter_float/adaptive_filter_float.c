/* ----------------------------------------------------
* Adaptive Filter FLOAT version
* ----------------------------------------------------
* 

*/

#include "adaptive_filter_float.h"
#include <limits.h>


//NLMS 1 STAGE WITH MULTI INPUT
void AdaptFilter_multi_init_float(AdaptFilter_multi_float *f)
{
    int i;
    int j;
    for (i = 0; i < ORDER; ++i){
        for(j=0; j < N_INPUT; j++){
            f->x[j][i] = 0.0f;
            f->weights[j][i] = 0.0f;
        }
    }
    f->d = 0.0f;

}

void AdaptFilter_multi_put_float(AdaptFilter_multi_float *f, float x_i0, float x_i1, float x_i2, float d_i)
{

    //With this I'm doing the shifting manually.
    //Potrei evitare questo: per rendere il tutto più efficiente
    f->d = d_i;
    for (int i = ORDER - 1; i > 0; i--) {
        for (int j=0; j<N_INPUT; j++){
            f->x[j][i] = f->x[j][i - 1];
        }
        }
    f->x[0][0] = x_i0;
    f->x[1][0] = x_i1;
    f->x[2][0] = x_i2;
    
}

float AdaptFilter_multi_get_float(AdaptFilter_multi_float *f, float mu)
{       

    float acc = 0.0f;
    for (int i = 0; i < ORDER; i++) {
        for(int j=0; j < N_INPUT; j++){
            acc += f->weights[j][i] * f->x[j][i];  // f->weights deve essere double[]
        }
    }

    float e = f->d - acc;  // f->d deve essere double

    float epsilon = 1e-8f;

    float norm = epsilon;
    for (int i = 0; i < ORDER; ++i) {
        for(int j=0; j < N_INPUT; j++){
            norm += f->x[j][i] * f->x[j][i];
        }
    }


    // Aggiorna i pesi usando NLMS
    for (int i = 0; i < ORDER; ++i) {
        for(int j=0; j < N_INPUT; j++){
            f->weights[j][i] += (mu/norm) * e * f->x[j][i];
        }
    }

    return e;
}


