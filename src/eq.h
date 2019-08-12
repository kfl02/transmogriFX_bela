#pragma once

#ifndef TRANSMOGRIFX_EQ_H
#define TRANSMOGRIFX_EQ_H

#include "constants.h"

enum eq_mode {
    PK_EQ,
    LOW_SHELF,
    HIGH_SHELF
};

typedef struct cx_t {
    float r;  //magnitude or //real
    float i;  //phase or //imaginary
} cx;

typedef struct eq_t {
    float b0;
    float b1;
    float b2;
    float a1;
    float a2;

    float y1;
    float y2;
    float x1;
    float x2;

    //Current settings
    eq_mode type;
    float fs;
    float f0;
    float Q;
    float G;

    //helper variables
    float c, s, alpha, w0, A;

} eq_coeffs;

typedef struct equalizer_t {
    size_t nbands;
    eq_coeffs **band;
} eq_filters;

eq_filters *make_equalizer(size_t nbands, float fstart_, float fstop_, float sample_rate);

void eq_compute_coeffs(eq_coeffs *cf, int type, float fs, float f0, float Q, float G);

void eq_update_gain(eq_coeffs *cf, float G);

float geq_tick(eq_filters *eq, float x_);

void geq_tick_n(eq_filters *eq, float *xn, size_t N);

void plot_response(float, float, int, eq_coeffs *, float, cx *);
//void plot_response(float f1, float f2, int pts, eq_coeffs* cf, float fs, cx *r)

#endif //TRANSMOGRIFX_EQ_H
