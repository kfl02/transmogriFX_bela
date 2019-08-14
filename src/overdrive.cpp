#include <stdlib.h>

#include <cmath>

#include "overdrive.h"

// Allocate the overdrive struct and set default values
overdrive *make_overdrive(unsigned int oversample, unsigned int bsz, float fs) {
    overdrive *od;

    od = (overdrive *) malloc(sizeof(overdrive));
    od->procbuf = (float *) malloc(sizeof(float) * bsz * oversample);

    for (int i = 0; i < bsz; i++) {
        od->procbuf[i] = 0.0f;
    }

    od->xn1 = 0.0f;
    od->xc1 = 0.0f;

    od->blksz = bsz;
    od->oversample = oversample;
    od->fs = fs;
    od->clipper_fs = ((float) oversample) * fs;
    od->inverse_oversample_float = 1.0f / ((float) oversample);

    // Set defaults
    od->gain = 30.0f;
    od->tone = 0.5f;
    od->level = 0.5f;
    od->dry = 1.0f;
    od->bypass = true;

    // Setup EQ stages
    compute_filter_coeffs_1p(&(od->anti_alias), LPF1P, od->clipper_fs,
                             4410.0f);  // 1/10th fs assuming 44.1k rate -- increasing filter order will help
    // down-play any aliasing artefacts, but increases CPU usage
    compute_filter_coeffs_1p(&(od->pre_emph), HPF1P, od->fs, 720.0f);
    compute_filter_coeffs_1p(&(od->post_emph), LPF1P, od->clipper_fs, 860.0f);
    compute_filter_coeffs_1p(&(od->tone_lp), LPF1P, od->fs, 1200.0f);
    compute_filter_coeffs_1p(&(od->tone_hp), HPF1P, od->fs, 1700.0f);

    return od;
}

void overdrive_cleanup(overdrive *od) {
    free(od->procbuf);
    free(od);
}

inline float sqr(float x) {
    return x * x;
}

//
// Clipping functions
//

// Quadratic clipping function
// Linear between nthrs and thrs, uses x - a*x^2 type of function above threshold
const float thrs = 0.8f;
const float nthrs = -0.72f;
const float f = 1.25f;

void clipper_tick(overdrive *od, int N, float *x, float *clean)  // Add in gain processing and dry mix
{
    float xn = 0.0f;
    float dx;
    float dc;
    float delta = 0.0f;
    float deltac = 0.0f;

    for (int i = 0; i < N; i++) {
        // Compute deltas for linear interpolation (upsampling)
        dx = (x[i] - od->xn1) * od->inverse_oversample_float;
        dc = (clean[i] - od->xc1) * od->inverse_oversample_float;
        od->xn1 = x[i];
        od->xc1 = clean[i];

        // Run clipping function at higher sample rate
        for (int n = 0; n < od->oversample; n++) {
            xn = x[i] + delta; // Linear interpolation up-sampling
            xn *= od->gain; // Apply gain
            clean[i] = clean[i] + deltac; // Upsample clean signal for mix
            delta += dx;
            deltac += dc;
            //Hard limiting
            if (xn >= 1.2f) {
                xn = 1.2f;
            }
            if (xn <= -1.12f) {
                xn = -1.12f;
            }

            //Soft clipping
            if (xn > thrs) {
                xn -= f * sqr(xn - thrs);
            }
            if (xn < nthrs) {
                xn += f * sqr(xn - nthrs);
            }

            // Pre-filter for down-sampling
            // Run de-emphasis and anti-aliasing filters
            xn = tick_filter_1p(&(od->post_emph), (od->dry * clean[i] + 0.7f * xn));
            xn = tick_filter_1p(&(od->anti_alias), xn);
        }

        // Reset linear interpolator state variables
        delta = 0.0f;
        deltac = 0.0f;

        // Zero-order hold downsampling assumes de-emphasis filter and anti-aliasing
        // filters sufficiently rejected harmonics > 1/2 base sample rate
        x[i] = xn;
    }
}

// Cubic clipping function
void cubic_clip(overdrive *od, int N, float asym, float *x, float *clean) {
    float xn;
    float dx = 0.0f;
    float dc = 0.0f;
    float delta = 0.0f;
    float deltac = 0.0f;

    for (unsigned int i = 0; i < N; i++) {
        // Compute deltas for linear interpolation (upsampling)
        xn = x[i] + delta; // Linear interpolation up-sampling
        xn *= od->gain; // Apply gain
        clean[i] = clean[i] + deltac; // Upsample clean signal for mix

        // TODO: never used
        delta += dx;
        deltac += dc;

        // Run clipping function at higher sample rate
        for (int n = 0; n < od->oversample; n++) {
            // Cubic clipping
            xn = xn * od->gain * 0.33f +
                 asym;  // Gain reduced because d/dx(x^3) = 3x ==> Small-signal gain of 3 built into the function

            if (xn <= -1.0f) {
                xn = -2.0f / 3.0f;
            } else if (xn >= 1.0f) {
                xn = 2.0f / 3.0f;
            } else {
                xn = xn - (1.0f / 3.0f) * xn * xn * xn;
            }

            // Pre-filter for down-sampling
            // Run de-emphasis and anti-aliasing filters
            xn = tick_filter_1p(&(od->post_emph), (od->dry * clean[i] + 0.7f * xn));
            xn = tick_filter_1p(&(od->anti_alias), xn);
        }
        // Reset linear interpolator state variables
        delta = 0.0f;
        deltac = 0.0f;

        // Zero-order hold downsampling assumes de-emphasis filter and anti-aliasing
        // filters sufficiently rejected harmonics > 1/2 base sample rate
        x[i] = xn;
    }
}

// Set EQ parameters to non-default
// Could be real-time user-configurable, but meant for
// configuring the type of overdrive
void od_set_cut_pre_emp(overdrive *od, float fc) {
    compute_filter_coeffs_1p(&(od->pre_emph), HPF1P, od->fs, fc);
}

void od_set_cut_post_emp(overdrive *od, float fc) {
    compute_filter_coeffs_1p(&(od->post_emph), HPF1P, od->fs, fc);
}

void od_set_cut_tone_lp(overdrive *od, float fc) {
    compute_filter_coeffs_1p(&(od->tone_lp), HPF1P, od->fs, fc);
}

void od_set_cut_tone_hp(overdrive *od, float fc) {
    compute_filter_coeffs_1p(&(od->tone_hp), HPF1P, od->fs, fc);
}

// Typical real-time user-configurable parameters
void od_set_drive(overdrive *od, float drive_db)   // 0 dB to 45 dB
{
    float drv = drive_db;

    if (drv < 0.0f) {
        drv = 0.0f;
    } else if (drv > 45.0f) {
        drv = 45.0f;
    }

    od->gain = std::pow(10.0f, drv / 20.0f);
}

void od_set_tone(overdrive *od, float hf_level_db) // high pass boost/cut, +/- 12dB
{
    float tone = hf_level_db;

    if (tone < -12.0f) {
        tone = -12.0f;
    } else if (tone > 12.0f) {
        tone = 12.0f;
    }

    od->tone = std::pow(10.0f, tone / 20.0f);
}

void od_set_level(overdrive *od, float outlevel_db) // -40 dB to +0 dB
{
    float vol = outlevel_db;

    if (vol < -40.0f) {
        vol = 40.0f;
    }
    if (vol > 0.0f) {
        vol = 0.0f;
    }

    od->level = std::pow(10.0f, vol / 20.0f);
}

void od_set_dry(overdrive *od, float dry) {
    if (dry < 0.0f) {
        od->dry = 0.0f;
    } else if (dry > 1.0f) {
        od->dry = 1.0f;
    } else {
        od->dry = dry;
    }
}

bool od_set_bypass(overdrive *od, bool bypass) {
    if (!bypass) {
        od->bypass = !od->bypass;
    } else {
        od->bypass = true;
    }

    return od->bypass;
}

// Run the overdrive effect
void overdrive_tick(overdrive *od, float *x) {
    const unsigned int n = od->blksz;

    if (od->bypass) {
        return;
    }

    // Run pre-emphasis filter
    for (unsigned int i = 0; i < n; i++) {
        od->procbuf[i] = tick_filter_1p(&(od->pre_emph), x[i]);
    }

    // Run the clipper
    clipper_tick(od, n, od->procbuf, x);  // Quadratic function: x - a*x^2
    //cubic_clip(od, n, 0.1, od->procbuf, x);   // Cubic: x - b*x^3

    // Output level and tone control
    for (unsigned int i = 0; i < n; i++) {
        x[i] = od->procbuf[i] * od->level;
        x[i] = 0.25f * x[i] + tick_filter_1p(&(od->tone_lp), x[i]) + od->tone * tick_filter_1p(&(od->tone_hp), x[i]);
    }
}
