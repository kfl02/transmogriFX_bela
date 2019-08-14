#include <cmath>

#include "iir_1pole.h"

//
// First order high-pass and low-pass filters.
//

void compute_filter_coeffs_1p(iir_1p *cf, iir_lp_mode type, float fs, float f0) {
    const float w0 = 2.0f * PI * f0 / fs;

    float a1;
    float b0, b1;
    float g = 1.0f;  // This could be brought out into a user-configurable param

    cf->gain = g;
    cf->fs = fs;

    switch (type) {
        default:
        case LPF1P:
            //1-pole high pass filter coefficients
            // H(z) = g * (1 - z^-1)/(1 - a1*z^-1)
            // Direct Form 1:
            //    h[n] = g * ( b0*x[n] - b1*x[n-1] ) - a1*y[n-1]
            // In below implementation gain is redistributed to the numerator:
            //    h[n] = gb0*x[n] - gb1*x[n-1] - a1*y[n-1]
            a1 = -expf(-w0);
            g = (1.0f + a1) / 1.12f;
            b0 = g;
            b1 = 0.12f * g; //0.12 zero improves RC filter emulation at higher freqs.
            break;

        case HPF1P:
            //1-pole high pass filter coefficients
            // H(z) = g * (1 - z^-1)/(1 - a1*z^-1)
            // Direct Form 1:
            //    h[n] = g * ( b0*x[n] -  b1*x[n-1] ) - a1*y[n-1]
            // In below implementation gain is redistributed to the numerator:
            //    h[n] = g*x[n] - g*x[n-1] - a1*y[n-1]
            a1 = -expf(-w0);
            g = (1.0f - a1) * 0.5f;
            b0 = g;
            b1 = -g;
            break;
    }

    cf->b0 = b0;
    cf->b1 = b1;
    cf->a1 = -a1;  // filter implementation uses addition instead of subtraction

    cf->y1 = 0.0;
    cf->x1 = 0.0;

    // Unused for a first order filter
    cf->a2 = 0.0f;
    cf->b2 = 0.0f;
    cf->x2 = 0.0f;
    cf->y2 = 0.0f;

}

void iir_get_response(iir_1p *cf, float n, float fstart, float fstop, float *frq, float *mag, float *phase) {
    const int RE = 0;
    const int IM = 1;

    // TODO: n is float already
    const float fn = n;
    const float fstep = (fstop - fstart) / fn;
    const float fs = cf->fs;

    float fi = fstart;
    float num[2] = {0, 0};
    float den[2] = {0, 0};

    // TODO: n is float
    for (int i = 0; i < n; i++) {
        num[RE] = cf->b0 + cf->b1 * cosf(2.0f * PI * fi / fs) + cf->b2 * cosf(4.0f * PI * fi / fs);
        num[IM] = -cf->b1 * sinf(2.0f * PI * fi / fs) - cf->b2 * sinf(4.0f * PI * fi / fs);
        den[RE] = 1.0f - cf->a1 * cosf(2.0f * PI * fi / fs) - cf->a2 * cosf(4.0f * PI * fi / fs);
        den[IM] = cf->a1 * sinf(2.0f * PI * fi / fs) + cf->a2 * sinf(4.0f * PI * fi / fs);

        mag[i] = 20.0f * std::log10(cf->gain *
                                       std::sqrt((num[RE] * num[RE] + num[IM] * num[IM]) /
                                     (den[RE] * den[RE] + den[IM] * den[IM])));
        phase[i] = 0.0f; // TODO: Implement phase calculation

        frq[i] = fi;
        fi += fstep;
    }
}


//
// Input s-domain biquad coefficients + filter gain
// Replaces num/den arrays with z-domain biquad coefficients
// Computed via Bilinear Transform
//
// sgain : Gain of s-domain transfer function or system linear gains to be 
//         distributed into z-domain coefficients
// fs_   : DSP system sampling rate
// kz_   : Frequency warping coefficient, or 0.0 for default
// num   : Biquad numerator:    num[2]*s^2 + num[1]*s + num[0]
// den   : Buquad denominator:  den[2]*s^2 + den[1]*s + den[0]
//
// den[1] and den[2] are negated for difference equation implementation in 
// MAC hardware (all operations summation).
//
// den[0] is always set to 0 for the z-domain biquad form.
//
// The intended difference equation implementation is of the following format:
// y[n] =   x[n-2]*num[2] + x[n-1]*num[1] + x[n]*num[0]
//        + y[n-2]*den[2] + y[n-1]*den[1]
//
// Default usage example:
//   s_biquad_to_z_biquad(sgain, fs, 0.0, num, den);
//

void s_biquad_to_z_biquad(float sgain, float fs_, float kz_, float *num, float *den) {
    float kz = kz_;
    const float fs = fs_;

    if (kz == 0.0f) {
        kz = 2.0f * fs;
    }

    const float kz2 = kz * kz;

    const float B0 = num[2] * sgain;
    const float B1 = num[1] * sgain;
    const float B2 = num[0] * sgain;
    const float A0 = den[2];
    const float A1 = den[1];
    const float A2 = den[0];

    num[2] = (B0 * kz2 + B1 * kz + B2) / (A0 * kz2 + A1 * kz + A2);
    num[1] = (2.0f * B2 - 2.0f * B0 * kz2) / (A0 * kz2 + A1 * kz + A2);
    num[0] = (B0 * kz2 - B1 * kz + B2) / (A0 * kz2 + A1 * kz + A2);

    den[2] = -(A0 * kz2 - A1 * kz + A2) / (A0 * kz2 + A1 * kz + A2);      // negated
    den[1] = -(2.0f * A2 - 2.0f * A0 * kz2) / (A0 * kz2 + A1 * kz + A2);    // negated
    den[0] = 0.0;  // should not be used
}
