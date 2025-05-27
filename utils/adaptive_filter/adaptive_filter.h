#ifndef ADAPTIVEFILTER_H_
#define ADAPTIVEFILTER_H_

#include "../../types.h"

#define ORDER 50


typedef struct AdaptFilter
{
  double x[ORDER]; //Buffer containing acc
  double weights[ORDER];
  double d;
  unsigned int last_index;
} AdaptFilter;

void AdaptFilter_init(AdaptFilter *f);
void AdaptFilter_put(AdaptFilter *f, double x, double d);
double AdaptFilter_get(AdaptFilter *f);

#endif

