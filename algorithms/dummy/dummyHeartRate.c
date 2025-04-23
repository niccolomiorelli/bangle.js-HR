/* ----------------------------------------------------
* DUMMY ALGORITHM
* -----------------------------------------------------
* Description:
* The algorithm return a fixed number (71) as output
*
*/

#include "../../types.h"
#include "dummyHeartRate.h"

/// Initialise step counting
void dummy_heartrate_init()
{
}

// process sample
int dummy_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz)
{
    return 710;
}
