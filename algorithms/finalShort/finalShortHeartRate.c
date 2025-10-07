/* ----------------------------------------------------
* FINAL SHORT ALGORITHM (simplified)
* ----------------------------------------------------
* This algorithm keeps only:
* - Band-pass filtering (FIR 0.6-3.5 Hz)
* - Standardization (z-score)
* - Spectral analysis with FFT (with zero padding) and peak selection
* Skips: NLMS adaptive filter and Kalman-based post-processing.
*/

#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>

#include "../../utils/bandpass_filter/bandpass_filter.h"
#include "../../utils/rolling_stats_float/rolling_stats_float.h"
#include "../../utils/fft_library_float/fft_library_float.h"
#include "../../types.h"

#include "finalShortHeartRate.h"

// --- PARAMETERS ---
#define WINDOW_STEP 64
#define SAMPLING_FREQ 25

// Physiological range (36-210 bpm)
#define LOW_FREQ_PHY_I 25
#define HIGH_FREQ_PHY_I 143

// --- STATE ---
static BPFilter bpFilter_ppg;
static Stats_float stats_ppg;

static float signal_PPG_buffer[WINDOW_LEN] = {0.0f};
static int signal_buffer_next_i = 0;
static int samples_since_last_HR = 0;

static float windowed_signal_in[WINDOW_LEN];
static float windowed_signal_out[WINDOW_LEN];
static complex_number_float fft_input[WINDOW_LEN];
static complex_number_float fft_input_padded[N_PAD];

// HR result (x10 bpm like the other algorithms)
static int HR = 0;
static bool first_window = true;

// simple circular buffer indexing helper
static inline int circ_idx(int head, int plus, int max) { return (head + plus) % max; }

// Moving average buffer for HR results (bpm*10)
#define HRM_AVG_LEN 5
static int hr_hist[HRM_AVG_LEN] = {0};
static int hr_hist_index = 0;
static int hr_hist_count = 0;

// --- INIT ---
void finalShort_heartrate_init()
{
    BPFilter_init(&bpFilter_ppg);
    rolling_stats_reset_float(&stats_ppg);

    for (int i = 0; i < WINDOW_LEN; ++i) {
        signal_PPG_buffer[i] = 0.0f;
        windowed_signal_in[i] = 0.0f;
        windowed_signal_out[i] = 0.0f;
    }
    signal_buffer_next_i = 0;
    samples_since_last_HR = 0;
    HR = 0;
    first_window = true;

    for (int i = 0; i < HRM_AVG_LEN; ++i) hr_hist[i] = 0;
    hr_hist_index = 0;
    hr_hist_count = 0;
}

// --- MAIN STEP ---
static int compute_fft_and_peak(void)
{
    // Prepare window from circular buffer
    for (int i = 0; i < WINDOW_LEN; i++) {
        int buffer_i = circ_idx(signal_buffer_next_i, i, WINDOW_LEN);
        windowed_signal_in[i] = signal_PPG_buffer[buffer_i];
    }

    // Hann window
    apply_hann_window_fftLib_float(windowed_signal_in, windowed_signal_out, WINDOW_LEN);

    // Copy to complex buffer
    for (int i = 0; i < WINDOW_LEN; i++) {
        fft_input[i].real = windowed_signal_out[i];
        fft_input[i].imag = 0.0f;
    }

    // Zero pad and FFT
    zero_pad_fftLib_float(fft_input, fft_input_padded, WINDOW_LEN);
    FFT_fftLib_float(fft_input_padded, 1.0f);

    // Peak selection in physiological band
    float max_fft_magnitude = 0.0f;
    int dominant_freq_index = 0;

    for (int i = LOW_FREQ_PHY_I; i <= HIGH_FREQ_PHY_I; i++) {
        float mag = (fft_input_padded[i].real * fft_input_padded[i].real +
                     fft_input_padded[i].imag * fft_input_padded[i].imag) / (float)N_PAD;
        // simple local-maximum check when possible
        if (i > LOW_FREQ_PHY_I && i < HIGH_FREQ_PHY_I) {
            float prev = (fft_input_padded[i-1].real * fft_input_padded[i-1].real +
                          fft_input_padded[i-1].imag * fft_input_padded[i-1].imag) / (float)N_PAD;
            float next = (fft_input_padded[i+1].real * fft_input_padded[i+1].real +
                          fft_input_padded[i+1].imag * fft_input_padded[i+1].imag) / (float)N_PAD;
            if (mag > prev && mag > next && mag > max_fft_magnitude) {
                max_fft_magnitude = mag;
                dominant_freq_index = i;
            }
        } else {
            if (mag > max_fft_magnitude) {
                max_fft_magnitude = mag;
                dominant_freq_index = i;
            }
        }
    }

    if (dominant_freq_index <= 0) {
        return HR; // keep last value if no peak found
    }

    // Convert index to bpm*10
    float freq_hz = (float)dominant_freq_index * (float)SAMPLING_FREQ / (float)N_PAD;
    int hr10 = (int)(freq_hz * 600.0f); // bpm*10

    // Update moving average buffer and compute average
    hr_hist[hr_hist_index] = hr10;
    hr_hist_index = (hr_hist_index + 1) % HRM_AVG_LEN;
    if (hr_hist_count < HRM_AVG_LEN) hr_hist_count++;

    long sum = 0;
    for (int i = 0; i < hr_hist_count; ++i) sum += hr_hist[i];
    int hr10_avg = (int)(sum / (hr_hist_count > 0 ? hr_hist_count : 1));
    return hr10_avg;
}

int finalShort_heartrate(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz)
{
    (void)accx; (void)accy; (void)accz; // unused in simplified version

    // Band-pass filter PPG
    BPFilter_put(&bpFilter_ppg, ppg);
    ppg_t ppg_filtered = BPFilter_get(&bpFilter_ppg);

    // Rolling z-score standardization
    rolling_stats_addValue_float((float)ppg_filtered, &stats_ppg);
    float mean_ppg = rolling_stats_get_mean_float(&stats_ppg);
    if (rolling_stats_get_variance_float(&stats_ppg) < 0.0f) {
        rolling_stats_reset_float(&stats_ppg);
        mean_ppg = 0.0f;
    }
    float std_ppg = rolling_stats_get_standard_deviation_float(&stats_ppg);
    if (std_ppg == 0.0f) std_ppg = 1.0f;
    float ppg_standardized = ((float)ppg_filtered - mean_ppg) / std_ppg;

    // Push into circular buffer
    signal_PPG_buffer[signal_buffer_next_i] = ppg_standardized;
    signal_buffer_next_i = (signal_buffer_next_i + 1) % WINDOW_LEN;
    samples_since_last_HR++;

    if (samples_since_last_HR == WINDOW_STEP) {
        if (first_window) {
            // first half window collected; wait for next
            first_window = false;
            samples_since_last_HR = 0;
        } else {
            HR = compute_fft_and_peak();
            samples_since_last_HR = 0;
        }
    }

    return HR;
}
