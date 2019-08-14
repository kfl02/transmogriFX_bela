#include <stdlib.h>

#include <cmath>

#include "biquad.h"

float *
make_butterworth_coeffs(int order) {
    float *coeffs = (float *) malloc(sizeof(float) * order);

    int n = order;
    float fn = (float) n;

    if (n % 2 == 0) {
        for (int k = 1; k <= n / 2; k++) {
            float fk = (float) k;

            coeffs[k - 1] =
                    1.0f /
                    (-2.0f * std::cos((2.0f * fk + fn - 1.0f) / (2.0f * fn) * PI)); //1/x returns filter stage Q factor
        }
    } else { //odd
        for (int k = 1; k <= (n - 1) / 2; k++) {
            float fk = (float) k;

            coeffs[k - 1] =
                    1.0f /
                    (-2.0f * std::cos((2.0f * fk + fn - 1.0f) / (2.0f * fn) * PI)); //1/x returns filter stage Q factor
        }
    }

    return coeffs;
}

void
biquad_update_coeffs(biquad_mode type, biquad_coeffs *cf, float fs, float f0, float Q) {
    const float w0 = 2.0f * PI * f0 / fs;
    const float c = std::cos(w0);
    const float s = std::sin(w0);
    const float alpha = s / (2.0f * Q);
    float a0, a1, a2;
    float b0, b1, b2;

    a0 = a1 = a2 = 0.0f;
    b0 = b1 = b2 = 0.0f;

    switch (type) {
        case LPF:
            b0 = (1.0f - c) / 2.0f;
            b1 = 1.0f - c;
            b2 = (1.0f - c) / 2.0f;
            a0 = 1.0f + alpha;
            a1 = -2.0f * c;
            a2 = 1.0f - alpha;
            break;

        case HPF:
            b0 = (1.0f + c) / 2.0f;
            b1 = -(1.0f + c);
            b2 = (1.0f + c) / 2.0f;
            a0 = 1.0f + alpha;
            a1 = -2.0f * c;
            a2 = 1.0f - alpha;
            break;

        case BPF:
            b0 = Q * alpha;
            b1 = 0.0f;
            b2 = -Q * alpha;
            a0 = 1.0f + alpha;
            a1 = -2.0f * c;
            a2 = 1.0f - alpha;
            break;

        default:
            break;
    }

    b0 /= a0;
    b1 /= a0;
    b2 /= a0;
    a1 /= a0;
    a2 /= a0;

    cf->b0 = b0;
    cf->b1 = b1;
    cf->b2 = b2;
    cf->a1 = a1;
    cf->a2 = a2;

    cf->Q = Q;
}

biquad_coeffs *
make_biquad(biquad_mode type, float fs, float f0, float Q) {
    biquad_coeffs *cf;

    cf = (biquad_coeffs *) malloc(sizeof(biquad_coeffs));

    biquad_update_coeffs(type, cf, fs, f0, Q);

    cf->y1 = 0.0f;
    cf->y2 = 0.0f;
    cf->x1 = 0.0f;
    cf->x2 = 0.0f;

    return cf;
}

float
run_filter(float x, biquad_coeffs *cf) {
    float y0 = cf->b0 * x + cf->b1 * cf->x1 + cf->b2 * cf->x2
               - cf->a1 * cf->y1 - cf->a2 * cf->y2;
    cf->x2 = cf->x1;
    cf->x1 = x;
    cf->y2 = cf->y1;
    cf->y1 = y0;

    return y0;
}
