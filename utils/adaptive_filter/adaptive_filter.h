#ifndef ADAPTIVEFILTER_H_
#define ADAPTIVEFILTER_H_

#include "../../types.h"

#define ORDER 75
#define N_INPUT 3

//FOR 1 STAGE SINGLE INPUT NLMS
typedef struct AdaptFilter1
{
  double x[ORDER]; //Buffer containing acc
  double weights[ORDER];
  double d;
  unsigned int last_index;
} AdaptFilter1;

void AdaptFilter1_init(AdaptFilter1 *f);
void AdaptFilter1_put(AdaptFilter1 *f, double x, double d);
double AdaptFilter1_get(AdaptFilter1 *f, double mu);

//FOR 1 STAGE MULTI INPUT NLMS
typedef struct AdaptFilter_multi
{
  double x[N_INPUT][ORDER]; //Buffer containing acc
  double weights[N_INPUT][ORDER];
  double d;
  unsigned int last_index;
} AdaptFilter_multi;

void AdaptFilter_multi_init(AdaptFilter_multi *f);
void AdaptFilter_multi_put(AdaptFilter_multi *f, double x_i0, double x_i1, double x_i2, double d_i);
double AdaptFilter_multi_get(AdaptFilter_multi *f, double mu);

#endif

