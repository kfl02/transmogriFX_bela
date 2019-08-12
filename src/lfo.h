#pragma once

#include "constants.h"

#ifndef TRANSMOGRIFX_LFO_H
#define TRANSMOGRIFX_LFO_H

enum lfo_mode {
    INT_TRI,
    TRI,
    SINE,
    SQUARE,
    EXP,
    RELAX,
    HYPER,
    HYPER_SINE
};

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
    lfo_mode lfo_type;
} lfoparams;

lfoparams *init_lfo(float, float, float);

void update_lfo(lfoparams *, float, float);

void set_lfo_type(lfoparams *, unsigned int);

void get_lfo_name(unsigned int, char *);

float run_lfo(lfoparams *);

#endif //TRANSMOGRIFX_LFO_H

