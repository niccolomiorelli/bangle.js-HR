
#include "rolling_stats.h"
#include <math.h>



void rolling_stats_reset(Stats *stats)
{
    stats->n = 0;
    stats->index = 0;
    stats->old_mean = 0.0;
    stats->mean = 0.0;
    stats->run_var = 0.0;
    stats->run_var_old = 0.0;
    for (int i=0;i<WINDOW_SIZE;i++){
        stats->window[i] = 0.0;
    }
    stats->full = 0;
}

void rolling_stats_addValue(double x, Stats *stats)
{
    int i = stats->index;
    double old_value = stats->window[i];

    stats->window[i] = x;
    stats->index = (i + 1) % WINDOW_SIZE;
    

    if (!stats->full)
    {
        stats->n += 1;
        double delta = x - stats->mean;
        stats->mean += delta / stats->n;
        stats->run_var += delta * (x - stats->mean);
        if (stats->n == WINDOW_SIZE){
            stats->full = 1;
        }
    }
    else
    {
        stats->old_mean = stats->mean;
        stats->mean += (x - old_value) / (double)WINDOW_SIZE;
        stats->run_var += (x + old_value - stats->old_mean - stats->mean)*(x - old_value);
        //I add this check:
        if (stats->run_var < 0.0){
            stats->run_var = stats->run_var_old;
        }
        stats->run_var_old = stats->run_var; 

    }
}


double rolling_stats_get_mean(Stats *stats)
{
    return stats->n ? (double)(stats->mean) : 0.0; //if stats->n != 0 return stats->mean else return 0.0

}

double rolling_stats_get_variance(Stats *stats)
{
    int denom = stats->full ? WINDOW_SIZE : stats->n; //if (stats->full) denom = WINDOW_SIZE else denom = stats->n
    return denom > 1 ? (double)(stats->run_var / ((double)(denom) - 1)) : 0.0; //if (denom > 1) return stats->run_var / (denom - 1) else return 0.0
}

double rolling_stats_get_standard_deviation(Stats *stats)
{
    double variance = rolling_stats_get_variance(stats);
    return (double)sqrt(variance);
}
