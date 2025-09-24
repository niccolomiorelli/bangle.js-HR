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
    double old_mean;
    double mean;
    double run_var;
    double run_var_old;
    double window[WINDOW_SIZE];
    int full;
} Stats;

void rolling_stats_reset(Stats *stats);

void rolling_stats_addValue(double x, Stats *stats);

double rolling_stats_get_mean(Stats *stats);

double rolling_stats_get_variance(Stats *stats);

double rolling_stats_get_standard_deviation(Stats *stats);

#endif