#ifndef ADAPTIVEFILTER_FLOAT_H_
#define ADAPTIVEFILTER_FLOAT_H_

#include "../../types.h"

#define ORDER 75
#define N_INPUT 3


//FOR 1 STAGE MULTI INPUT NLMS
typedef struct AdaptFilter_multi_float
{
  float x[N_INPUT][ORDER]; //Buffer containing acc
  float weights[N_INPUT][ORDER];
  float d;
  unsigned int last_index;
} AdaptFilter_multi_float;

void AdaptFilter_multi_init_float(AdaptFilter_multi_float *f);
void AdaptFilter_multi_put_float(AdaptFilter_multi_float *f, float x_i0, float x_i1, float x_i2, float d_i);
float AdaptFilter_multi_get_float(AdaptFilter_multi_float *f, float mu);

#endif