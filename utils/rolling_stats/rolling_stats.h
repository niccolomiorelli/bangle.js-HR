// compute statistics on a rolling basis.
// Based on the Welford algorithm, see https : // www.johndcook.com/blog/standard_deviation/

#ifndef ROLLING_STATS
#define ROLLING_STATS

#include <math.h>

#define WINDOW_SIZE 128

typedef struct Stats
{
    unsigned long n;
    int index;
    float old_mean;
    float mean;
    float run_var;
    float window[WINDOW_SIZE];
    int full;
} Stats;

void rolling_stats_reset(Stats *stats);

void rolling_stats_addValue(float x, Stats *stats);

float rolling_stats_get_mean(Stats *stats);

float rolling_stats_get_variance(Stats *stats);

float rolling_stats_get_standard_deviation(Stats *stats);

#endif