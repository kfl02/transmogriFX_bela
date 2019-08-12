#pragma once

#ifndef TRANSMOGRIFX_BIQUAD_H
#define TRANSMOGRIFX_BIQUAD_H

#include "constants.h"

enum biquad_mode {
    LPF,
    HPF,
    BPF
};

typedef struct biquad_t {
    float b0;
    float b1;
    float b2;
    float a1;
    float a2;

    float y1;
    float y2;
    float x1;
    float x2;

    //operational parameters
    float Q;
} biquad_coeffs;


biquad_coeffs *make_biquad(biquad_mode, float, float, float);
//make_biquad(int type, biquad_coeffs* cf, float fs, float f0, float Q);

void biquad_update_coeffs(biquad_mode , biquad_coeffs *, float, float, float);
// biquad_update_coeffs(int type, biquad_coeffs* cf, float fs, float f0, float Q);

//get an array of Q factor values for an N order butterworth filter.
float *make_butterworth_coeffs(int);

float run_filter(float, biquad_coeffs *);
//run_filter(float x, biquad_coeffs* cf);

#endif //TRANSMOGRIFX_BIQUAD_H
