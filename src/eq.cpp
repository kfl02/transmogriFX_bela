#include <stdlib.h>

#include <cmath>

#include "eq.h"

void eq_compute_coeffs(eq_coeffs *cf, eq_mode type, float fs, float f0, float Q, float G) {
    const float A = std::pow(10.0f, G / 40.0f);
    const float w0 = 2.0f * PI * f0 / fs;
    const float c = std::cos(w0);
    const float s = std::sin(w0);
    const float alpha = s / (2.0f * Q);

    float a0, a1, a2;
    float b0, b1, b2;

    a0 = a1 = a2 = 0.0f;
    b0 = b1 = b2 = 0.0f;

    cf->fs = fs;
    cf->f0 = f0;
    cf->Q = Q;
    cf->G = G;
    cf->type = type;

    cf->A = A;
    cf->w0 = w0;
    cf->c = c;
    cf->s = s;
    cf->alpha = alpha;

    switch (type) {
        case PK_EQ:
            b0 = 1.0f + alpha * A;
            b1 = -2.0f * c;
            b2 = 1.0f - alpha * A;
            a0 = 1.0f + alpha / A;
            a1 = -2.0f * c;
            a2 = 1.0f - alpha / A;
            break;

        case LOW_SHELF:
            b0 = A * ((A + 1.0f) - (A - 1.0f) * c + 2.0f * std::sqrt(A) * alpha);
            b1 = 2.0f * A * ((A - 1.0f) - (A + 1.0f) * c);
            b2 = A * ((A + 1.0f) - (A - 1.0f) * c - 2.0f * std::sqrt(A) * alpha);
            a0 = (A + 1.0f) + (A - 1.0f) * c + 2.0f * std::sqrt(A) * alpha;
            a1 = -2.0f * ((A - 1.0f) + (A + 1.0f) * c);
            a2 = (A + 1.0f) + (A - 1.0f) * c - 2.0f * std::sqrt(A) * alpha;
            break;

        case HIGH_SHELF:
            b0 = A * ((A + 1.0f) + (A - 1.0f) * c + 2.0f * std::sqrt(A) * alpha);
            b1 = -2.0f * A * ((A - 1.0f) + (A + 1.0f) * c);
            b2 = A * ((A + 1.0f) + (A - 1.0f) * c - 2.0f * std::sqrt(A) * alpha);
            a0 = (A + 1.0f) - (A - 1.0f) * c + 2.0f * std::sqrt(A) * alpha;
            a1 = 2.0f * ((A - 1.0f) - (A + 1.0f) * c);
            a2 = (A + 1.0f) - (A - 1.0f) * c - 2.0f * std::sqrt(A) * alpha;
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
    //negate 'a' coefficients so addition can be used in filter
    cf->a1 = -a1;
    cf->a2 = -a2;

}

void eq_update_gain(eq_coeffs *cf, float G) {
    const float A = std::pow(10.0f, G / 40.0f);
    const float iA = 1.0f / A;
    const float sqrtA2 = 2.0f * std::sqrt(A);

    float a0, a1, a2;
    float b0, b1, b2;

    a0 = a1 = a2 = 0.0f;
    b0 = b1 = b2 = 0.0f;

    const float c = cf->c;
    const float alpha = cf->alpha;

    switch (cf->type) {
        case PK_EQ:
            b0 = 1.0f + alpha * A;
            b1 = -2.0f * c;
            b2 = 1.0f - alpha * A;
            a0 = 1.0f + alpha * iA;
            a1 = -2.0f * c;
            a2 = 1.0f - alpha * iA;
            break;

        case LOW_SHELF:
            b0 = A * ((A + 1.0f) - (A - 1.0f) * c + sqrtA2 * alpha);
            b1 = 2.0f * A * ((A - 1.0f) - (A + 1.0f) * c);
            b2 = A * ((A + 1.0f) - (A - 1.0f) * c - sqrtA2 * alpha);
            a0 = (A + 1.0f) + (A - 1.0f) * c + sqrtA2 * alpha;
            a1 = -2.0f * ((A - 1.0f) + (A + 1.0f) * c);
            a2 = (A + 1.0f) + (A - 1.0f) * c - sqrtA2 * alpha;
            break;

        case HIGH_SHELF:
            b0 = A * ((A + 1.0f) + (A - 1.0f) * c + sqrtA2 * alpha);
            b1 = -2.0f * A * ((A - 1.0f) + (A + 1.0f) * c);
            b2 = A * ((A + 1.0f) + (A - 1.0f) * c - sqrtA2 * alpha);
            a0 = (A + 1.0f) - (A - 1.0f) * c + sqrtA2 * alpha;
            a1 = 2.0f * ((A - 1.0f) - (A + 1.0f) * c);
            a2 = (A + 1.0f) - (A - 1.0f) * c - sqrtA2 * alpha;
            break;

        default:
            break;
    }

    float ia0 = 1.0f / a0;

    b0 *= ia0;
    b1 *= ia0;
    b2 *= ia0;
    a1 *= ia0;
    a2 *= ia0;

    cf->b0 = b0;
    cf->b1 = b1;
    cf->b2 = b2;
    //negate 'a' coefficients so addition can be used in filter
    cf->a1 = -a1;
    cf->a2 = -a2;
}

eq_coeffs *make_eq_band(eq_mode type, eq_coeffs *cf, float fs, float f0, float Q, float G) {
    cf = (eq_coeffs *) malloc(sizeof(eq_coeffs));

    eq_compute_coeffs(cf, type, fs, f0, Q, G);

    cf->y1 = 0.0;
    cf->y2 = 0.0;
    cf->x1 = 0.0;
    cf->x2 = 0.0;

    return cf;
}

eq_filters *make_equalizer(size_t nbands, float fstart_, float fstop_, float sample_rate) {
    eq_filters *eq;

    const float G = 0.0f;

    float f1 = fstop_;
    float f0 = fstart_;

    if (f1 > sample_rate / 2.0f) {
        f1 = sample_rate / 2.0f;
    }

    const float nb = (float) nbands;
    const float m = std::pow(2.0f, ((logf(f1 / f0) / logf(2.0f)) / (nb - 1.0f)));
    const float Q = 0.25f * (m + 1.0f) / (m - 1.0f);

    eq = (eq_filters *) malloc(sizeof(eq_filters));
    eq->band = (eq_coeffs **) malloc((nbands + 2) * sizeof(eq_coeffs *));
    eq->nbands = nbands;

    for (int i = 0; i < nbands; i++) {
        //printf("f0 = %f\tQ = %f\tfl=%f\tfh=%f\n", f0, Q, f0-0.5*f0/Q, f0+0.5*f0/Q);
        eq->band[i] = make_eq_band(PK_EQ, eq->band[i], sample_rate, f0, Q, G);
        eq->band[i]->f0 = f0;
        eq->band[i]->Q = Q;
        eq->band[i]->G = G;
        eq->band[i]->fs = sample_rate;
        eq->band[i]->type = PK_EQ;

        f0 *= m;
    }
    //Create low shelf and high shelf filters
    eq->band[nbands] = make_eq_band(LOW_SHELF, eq->band[nbands], sample_rate, 336.0f, 0.5f, G);
    eq->band[nbands + 1] = make_eq_band(HIGH_SHELF, eq->band[nbands + 1], sample_rate, 1418.0f, 0.5f, G);

    eq->band[nbands]->f0 = 336.0f;
    eq->band[nbands]->Q = Q;
    eq->band[nbands]->G = G;
    eq->band[nbands]->fs = sample_rate;
    eq->band[nbands]->type = LOW_SHELF;

    eq->band[nbands + 1]->f0 = 1418.0f;
    eq->band[nbands + 1]->Q = Q;
    eq->band[nbands + 1]->G = G;
    eq->band[nbands + 1]->fs = sample_rate;
    eq->band[nbands + 1]->type = HIGH_SHELF;

    return eq;
}

inline float tick_eq_band(eq_coeffs *cf, float x) {
    float y0 = cf->b0 * x + cf->b1 * cf->x1 + cf->b2 * cf->x2
               + cf->a1 * cf->y1 + cf->a2 * cf->y2;
    cf->x2 = cf->x1;
    cf->x1 = x;
    cf->y2 = cf->y1;
    cf->y1 = y0;

    return y0;
}


float geq_tick(eq_filters *eq, float x_) {
    float x = x_;

    for (size_t i = 0; i < (eq->nbands + 2); i++) {
        x = tick_eq_band(eq->band[i], x);
    }

    return x;
}

void geq_tick_n(eq_filters *eq, float *xn, size_t N) {
    const size_t cnt = eq->nbands + 2;
    float x;

    for (size_t k = 0; k < N; k++) {
        x = xn[k];

        for (size_t i = 0; i < cnt; i++) {
            x = tick_eq_band(eq->band[i], x);
        }

        xn[k] = x;
    }
}

inline float sqr(float x) {
    return x * x;
}

void plot_response(float f1, float f2, int pts, eq_coeffs *cf, float fs, cx *r) {
    const float dp = 2.0f * PI;
    // z^-1 = cos(w) - j*sin(w)
    // z^-2 = cos(2*w) - j*sin(2w)

    float num[2], den[2];

    num[0] = num[1] = 0.0f;
    den[0] = den[1] = 0.0f;

    const float dw = (dp * (f2 - f1) / (float) pts) / fs;
    const int n = std::abs(pts) + 1;

    float w = dp * f1 / fs;

    for (int i = 0; i < n; i++) {
        float magn;
        float magd;
        float an;
        float ad;

        num[0] = cf->b0 + cf->b1 * std::cos(-w) + cf->b2 * std::cos(-2.0f * w); //numerator real part
        num[1] = cf->b1 * std::sin(-w) + cf->b2 * std::sin(-2.0f * w);
        den[0] = 1.0f + cf->a1 * std::cos(-w) + cf->a2 * std::cos(-2.0f * w); //numerator real part
        den[1] = cf->a1 * std::sin(-w) + cf->a2 * std::sin(-2.0f * w);

        magn = std::sqrt(sqr(num[0]) + sqr(num[1]));
        magd = std::sqrt(sqr(den[0]) + sqr(den[1]));
        an = atanf(num[1] / num[0]); //angle of numerator
        ad = atanf(den[1] / den[0]); //angle of denominator

        r[i].r = magn / magd;
        r[i].i = an - ad;

        //printf("%f\t%f\t%f\n", fs*w/dp, 20.0*log10(r[i].r), 180.0*r[i].i/(PI));

        w += dw;
    }
}
