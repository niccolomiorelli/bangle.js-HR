#ifndef ADAPTIVEFILTER_H_
#define ADAPTIVEFILTER_H_

#include "../../types.h"

#define ORDER 50


typedef struct AdaptFilter
{
  float x[ORDER]; //Buffer containing acc
  float weights[ORDER];
  float d;
  unsigned int last_index;
} AdaptFilter;

void AdaptFilter_init(AdaptFilter *f);
void AdaptFilter_put(AdaptFilter *f, float x, float d);
float AdaptFilter_get(AdaptFilter *f);

#endif

