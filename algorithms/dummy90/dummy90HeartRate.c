/* ----------------------------------------------------
* DUMMY90 ALGORITHM
* -----------------------------------------------------
* Description:
* The algorithm returns a fixed number (900) as output,
* representing 90 BPM heart rate
*
*/

#include "../../types.h"
#include "dummy90HeartRate.h"

/// Initialise heart rate algorithm
void dummy90_heartrate_init()
{
    // No initialization needed for this dummy algorithm
}

// process sample
int dummy90_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz)
{
    return 900;  // Always return 900 (90 BPM)
}