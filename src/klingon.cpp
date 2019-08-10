#include <math.h>
#include <stdlib.h>

#include "klingon.h"

// Allocate the klingon struct and set default values
klingon *make_klingon(klingon *kot, unsigned int oversample, unsigned int bsz, float fs) {
    kot = (klingon *) malloc(sizeof(klingon));
    kot->procbuf = (float *) malloc(sizeof(float) * bsz * oversample);

    // Set up interpolation table
    const char *soft_file = "soft_clip.txt";
    const char *hard_file = "hard_clip.txt";
    const char *limit_file = "output_limit.txt";

    load_vi_data(&(kot->clip), (char *) soft_file);
    load_vi_data(&(kot->hard_clip), (char *) hard_file);
    load_vi_data(&(kot->output_limit), (char *) limit_file);

    for (int i = 0; i < bsz; i++) {
        kot->procbuf[i] = 0.0f;
    }

    kot->xn1 = 0.0f;

    kot->blksz = bsz;
    kot->oversample = oversample;
    kot->fs = fs;
    kot->clipper_fs = ((float) oversample) * fs;
    kot->inverse_oversample_float = 1.0 / ((float) oversample);

    // Set defaults
    kot->gain = 50.0f;
    kot->tone = 0.5f; // Not used -- initialized in tonestack
    kot->hard = 0.0f;
    kot->level = 0.5f;
    kot->bypass = true;

    // Setup EQ stages
    compute_filter_coeffs_1p(&(kot->anti_alias), LPF1P, kot->clipper_fs,
                             kot->fs / 4.0);  // down-play any aliasing artefacts
    compute_filter_coeffs_1p(&(kot->pre_emph482), HPF1P, kot->fs, 482.29f);
    compute_filter_coeffs_1p(&(kot->pre_emph589), HPF1P, kot->fs, 589.46f);

    // First stage gains
    // ( g482*pre_emph482(x) + g589*pre_emph589(x) )*gain + x
    kot->g482 = 1.0f / 33.0f;  //1/33k resistor ratio
    kot->g589 = 1.0f / 27.0f;  //1/27k resistor ratio
    // 1rst stage gain will be 10 to 110

    // First-stage pre-emphasis, if using other schematic sources
    // (biquad_pre_emph)*gain + x
    float k = 1000.0f;
    float n = 1e-9f;
    float r1 = 27.0f * k;
    float r2 = 33.0f * k;
    float c1 = 10.0f * n;
    float c2 = 10.0f * n;

    float num[3];
    float den[3];

    compute_s_biquad(r1, r2, c1, c2, num, den);      // s-domain coefficients
    s_biquad_to_z_biquad(1.0f, fs, 0.0f, num, den);    // compute bilinear transform

    kot->pre_emph_biquad.a1 = den[1];
    kot->pre_emph_biquad.a2 = den[2];
    kot->pre_emph_biquad.b0 = num[0];
    kot->pre_emph_biquad.b1 = num[1];
    kot->pre_emph_biquad.b2 = num[2];
    kot->pre_emph_biquad.x1 = 0.0f;
    kot->pre_emph_biquad.x2 = 0.0f;
    kot->pre_emph_biquad.y1 = 0.0f;
    kot->pre_emph_biquad.y2 = 0.0f;
    kot->pre_emph_biquad.gain = 1.0f;
    kot->pre_emph_biquad.fs = kot->fs;

    // Setup KOT tonestack
    kotstack_init(&(kot->stack), kot->fs);

    // Next compute filter cut-offs and inter-stage gains as dependent upon
    // gain pot rotational postion
    //

    // Extract pot resistance from gain (dB) setting
    float pot_fb = (kot->gain - 10.0f);  // Portion in op amp feedback
    float pot_ff = 100.0f - pot_fb + 10.0f; // Portion in series with following gain stage input

    float c_fb = 100e-12f;  // Sets low-pass roll-off
    float c_ff = 100e-9f;   // Sets high-pass feeding into second gain stage

    float fc_fb = 1.0f / (2.0f * PI * 10000.0f * c_fb);  // Low pass cut-off
    float fc_ff = 1.0f / (2.0f * PI * 1000.0f * pot_ff * c_ff);     // High-pass cut-off

    compute_filter_coeffs_1p(&(kot->pre_emph159), HPF1P, kot->fs, fc_ff);
    compute_filter_coeffs_1p(&(kot->post_emph), LPF1P, kot->clipper_fs, fc_fb);

    //Second stage gains
    kot->g159 = 1.0e-3f / pot_ff;  // Current fed into second stage amp 1/(10k + (1-x)*100k), ratio of 100k gain pot

    return kot;
}

void klingon_cleanup(klingon *kot) {
    vi_trace_cleanup(&(kot->clip));
    free(kot->procbuf);
    free(kot);
}

inline float sqr(float x) {
    return x * x;
}


//
// King of Tone, Marshall Bluesbreaker, etc.
//  First pre-emphasis stage high-pass second-order filter response
//    (VIN-)-*--/\/\/\/---||---*---||---(Io->)--GND
//           |    R2      C2   |   C1
//           *-----/\/\/\/-----*
//                   R1
// Transfer function is Io(s)/Vin(s)
//
// Default usage example:
//   compute_s_biquad(r1, r2, c1, c2, num, den);
//

void compute_s_biquad(float r1, float r2, float c1, float c2, float *num, float *den) {
    float ga = (r1 + r2) / (r1 * r2);
    float gs = 1.0f / ga;
    float z0 = 1.0f / (r2 * c2);
    float p0 = 1.0f / ((r1 + r2) * c2);

    num[0] = 0.0f;
    num[1] = ga * p0;
    num[2] = ga;

    den[0] = p0 / (c1 * gs);
    den[1] = (c1 * gs * z0 + 1.0f) / (c1 * gs);
    den[2] = 1.0f;

}

//
// Clipping functions
//  Nonlinear function applied by Lagrange interpolation of look-up mutable
//  Voltage/Current function (usually exported from SPICE)
//


void clipper_tick(klingon *kot, int N, float *x, float *clean)  // Add in gain processing and dry mix
{
    float xn = 0.0f;
    float dx = 0.0f;
    float delta = 0.0f;
    float tmp = 0.0f;

    for (int i = 0; i < N; i++) {
        // Compute deltas for linear interpolation (upsampling)
        dx = (x[i] - kot->xn1) * kot->inverse_oversample_float;

        // Run clipping function at higher sample rate
        for (int n = 0; n < kot->oversample; n++) {
            xn = tick_filter_1p(&(kot->post_emph), kot->xn1 + delta); // Linear interpolation up-sampling
            delta += dx;

            if (xn > 300.0e-6f) {
                xn = 300.0e-6f;
            }
            if (xn < -300.0e-6f) {
                xn = -300.0e-6f;
            }

            tmp = xn;

            // Run nonlinear function defined from text file
            xn = vi_trace_interp(&(kot->clip), tmp);
            tmp = vi_trace_interp(&(kot->hard_clip), tmp);

            xn = kot->hard * tmp + (1.0f - kot->hard) * xn;

            // Run anti-aliasing filter
            xn = tick_filter_1p(&(kot->anti_alias), xn);
        }

        // Reset linear interpolator state variables
        kot->xn1 = x[i];
        delta = 0.0f;

        // Zero-order hold downsampling assumes de-emphasis filter and anti-aliasing
        // filters sufficiently rejected harmonics > 1/2 base sample rate
        x[i] = xn;
    }
}

// Typical real-time user-configurable parameters
void kot_set_drive(klingon *kot, float drive_db)   // 0 dB to 45 dB
{
    float drv = drive_db;

    if (drv < 0.0f) {
        drv = 0.0f;
    } else if (drv > 45.0f) {
        drv = 45.0f;
    }

    // Convert gain given in dB to absolute value
    drv = powf(10.0f, (drv + 0.341f) / 20.0f);

    // Work backward through gain stages to get pot value adjustment
    float gl = 33.0f * 27.0f / (33.0f + 27.0f);
    float gb = 22.0f;

    // pot setting required to get requested gain
    kot->gain = (drv / gb - 1.0f) * gl;

    //
    // Next compute filter cut-offs and inter-stage gains as dependent upon
    // gain pot rotational postion
    //

    // Extract pot resistance from gain (dB) setting
    float pot_fb = (kot->gain - 10.0f);  // Portion in op amp feedback
    float pot_ff = 100.0f - pot_fb + 10.0f; // Portion in series with following gain stage input
    float c_ff = 100e-9f;   // Sets high-pass feeding into second gain stage
    float fc_ff = 1.0f / (2.0f * PI * 1000.0f * pot_ff * c_ff);     // High-pass cut-off

    compute_filter_coeffs_1p(&(kot->pre_emph159), HPF1P, kot->fs, fc_ff);

    //Second stage gains
    kot->g159 = 1.0e-3f / pot_ff;  // 220k/(10k + (1-x)*100k), ratio of 100k gain pot
}

//
// Tone Knob setting
//   High-pass cut
// 
void kot_set_tone(klingon *kot, float hf_level_db) {
    // Range limiting and converter from logarithmic function
    // to linear function of pot resistance ratios used for
    // computing tonestack IIR coefficients
    float tone = hf_level_db;

    if (tone < -60.0f) {
        tone = -60.0f;
    } else if (tone > 0.0f) {
        tone = 0.0f;
    }

    kot->tone = powf(10.0f, tone / 20.0f);  // not used elsewhere, might delete

    // Tone control moved to kot_tonestack.h/cpp 
    kotstack_set_tone(&(kot->stack), kot->tone);
}

void kot_set_boost(klingon *kot, float boost) // second high pass cut
{
    // kotstack_set_boost handles out-of-range error checking
    kotstack_set_boost(&(kot->stack), boost);
}

void kot_set_mix(klingon *kot, float hard)  // Dry/Wet control, 0.0 to 1.0
{
    float mix = hard;

    if (hard > 1.0f) {
        mix = 1.0f;
    } else if (hard < 0.0f) {
        mix = 0.0f;
    } else {
        mix = hard;
    }

    kot->hard = mix;
}

void kot_set_level(klingon *kot, float outlevel_db) // -40 dB to +0 dB
{
    float vol = outlevel_db;

    if (vol < -40.0f) {
        vol = -40.0f;
    }
    if (vol > 0.0f) {
        vol = 0.0f;
    }

    kot->level = powf(10.0f, vol / 20.0f);
}

bool kot_set_bypass(klingon *kot, bool bypass) {
    if (!bypass) {
        kot->bypass = !kot->bypass;
    } else {
        kot->bypass = true;
    }

    return kot->bypass;
}

// Run the klingon effect
void klingon_tick(klingon *kot, float *x) {
    // TODO: resolve unsigned int vs. int
    unsigned int n = kot->blksz;

    if (kot->bypass) {
        return;
    }

    // Run pre-emphasis filters
    for (unsigned int i = 0; i < n; i++) {
        // "Matsumin" schematic rendering (likely incorrect)
        //kot->procbuf[i] = tick_filter_1p(&(kot->pre_emph589), kot->g589*x[i]);
        //kot->procbuf[i] += tick_filter_1p(&(kot->pre_emph482), kot->g482*x[i]);

        // Other internet source schematic rendering (likely correct due to similarity with Marshall Bluesbreaker)
        kot->procbuf[i] = tick_filter_biquad(&(kot->pre_emph_biquad), x[i]);

        // First stage gain
        kot->procbuf[i] *= kot->gain;
        kot->procbuf[i] += x[i];

        // Second stage pre-emphasis
        x[i] = tick_filter_1p(&(kot->pre_emph159), kot->procbuf[i]);
        kot->procbuf[i] = x[i] * kot->g159;  // Current fed into second op amp inverting terminal
    }

    // Run the clipper
    clipper_tick(kot, n, kot->procbuf, x);  //

    // Output level and tone control
    for (unsigned int i = 0; i < n; i++) {
        x[i] = kot->level * kotstack_tick(&(kot->stack), kot->procbuf[i]);
        x[i] = vi_trace_interp(&(kot->output_limit),
                               x[i]); //limit to -1.0 to 1.0, but more gracefully than digital clip/limit
    }
}
