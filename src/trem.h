#pragma once

#ifndef TRANSMOGRIFX_TREM_H
#define TRANSMOGRIFX_TREM_H

#include "lfo.h"

const int MODERATE = 0;
const int TREM_SQUARE = 1;

typedef struct trem_t {
    float fs;
    bool bypass;

    float depth;
    float gain;
    int lfo_type;

    //modulation
    lfoparams *lfo;

} trem_coeffs;

//Initialize the struct and allocate memory
trem_coeffs *make_trem(float);

// input sample, pot gain (0...1), iwah struct
float trem_tick(trem_coeffs *, float);

void trem_tick_n(trem_coeffs *, float *, int);

//Select preset circuit voicings
void trem_circuit_preset(trem_coeffs *, int);

//settings
void trem_set_lfo_rate(trem_coeffs *, float);

void trem_set_lfo_depth(trem_coeffs *, float);

void trem_set_lfo_gain(trem_coeffs *, float);

void trem_set_lfo_type(trem_coeffs *, unsigned int);

bool trem_toggle_bypass(trem_coeffs *);

#endif //TRANSMOGRIFX_TREM_H
