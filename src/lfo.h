#pragma once

#include "constants.h"

#ifndef TRANSMOGRIFX_LFO_H
#define TRANSMOGRIFX_LFO_H

const int INT_TRI = 0;
const int TRI = 1;
const int SINE = 2;
const int SQUARE = 3;
const int EXP = 4;
const int RELAX = 5;
const int HYPER = 6;
const int HYPER_SINE = 7;
const int MAX_LFOS = 7;


typedef struct lfoparams_t {
    //Integrated triangle wave state variables
    float k, k_;
    float nk, nk_;
    float psign, psign_;
    float nsign, nsign_;
    float sign, sign_;
    float lfo;
    float x;
    int startup_delay;

    //Triangle wave state variables
    float ktri;
    float trisign;
    float trilfo;

    //Sine oscillator state variables
    float ksin;
    float sin_part;
    float cos_part;

    //Relaxation oscillator parameters
    float rlx_sign;
    float rlx_min;
    float rlx_max;
    float rlx_k, rlx_ik;
    float rlx_lfo;

    //Exponential oscillator parameters
    float exp_k;
    float exp_ik;
    float exp_min;
    float exp_max;
    float exp_sv;
    float exp_x;

    //globals
    float current_rate;

    //Select LFO type 
    unsigned int lfo_type;
} lfoparams;

lfoparams *init_lfo(lfoparams *, float, float, float);

void update_lfo(lfoparams *, float, float);

void set_lfo_type(lfoparams *, unsigned int);

void get_lfo_name(unsigned int, char *);

float run_lfo(lfoparams *);

#endif //TRANSMOGRIFX_LFO_H

