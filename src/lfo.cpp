#include <stdio.h>
#include <stdlib.h>

#include <cmath>

#include "lfo.h"

lfoparams *init_lfo(float fosc, float fs, float phase) {
    lfoparams *lp;

    lp = (lfoparams *) malloc(sizeof(lfoparams));

    const float ts = 1.0f / fs;
    const float frq = 2.0f * fosc;
    const float t = 4.0f * frq * frq * ts * ts;

    float p = phase / 180.0f;  //Phase can be as large as desired to delay LFO startup

    if (p < 0.0f) {
        p *= -1.0f;
    }  //don't let it be negative

    p /= frq;
    p *= fs;

    const int phdly = (int) p;

    lp->startup_delay = phdly;

    //Integrated triangle wave LFO (quasi-sinusoidal, linear pitch shift)
    lp->k_ = lp->k = 2.0f * ts * frq;
    lp->nk_ = lp->nk = -2.0f * ts * frq;
    lp->psign_ = lp->psign = t;
    lp->nsign_ = lp->nsign = -t;
    lp->sign_ = lp->sign = t;
    lp->lfo = 0.0f;
    lp->x = 0.0f;

    //Triangle wave LFO variables
    lp->ktri = frq / fs;
    lp->trisign = 1.0f;
    p = frq * phase / (360.0f * fosc);

    if (p >= 1.0f) {
        p -= 1.0f;
        lp->trisign = -1.0f;
    }
    if (p < 0.0f) {
        p = 0.0f;
        lp->trisign = 1.0f;
    }

    lp->trilfo = p;

    //Sine wave LFO variables
    lp->ksin = PI * frq / fs;
    lp->sin_part = sin(2.0f * PI * phase / 360.0f);
    lp->cos_part = cos(2.0f * PI * phase / 360.0f);

    //Relaxation oscillator parameters
    static const float ie = 1.0f / (1.0f - 1.0f / E);
    const float k = expf(-2.0f * fosc / fs);

    lp->rlx_k = k;
    lp->rlx_ik = 1.0f - k;

    lp->rlx_sign = ie;
    lp->rlx_max = ie;
    lp->rlx_min = 1.0f - ie;

    lp->rlx_lfo = 0.0;

    //Exponential oscillator parameters
    k = expf(-2.0f * 1.3133f * fosc / fs);
    lp->exp_ik = k;
    lp->exp_k = 1.0f / k;
    lp->exp_x = k;
    lp->exp_min = 1.0f / E;
    lp->exp_max = 1.0f + 1.0f / E;
    lp->exp_sv = lp->exp_min;

    //Globals
    lp->current_rate = fosc;
    // TODO: type should be const/enum
    lp->lfo_type = 0; //integrated triangle

    return lp;
}

void update_lfo(lfoparams *lp, float fosc, float fs) {
    const float ts = 1.0f / fs;
    const float frq = 2.0f * fosc;
    const float t = 4.0f * frq * frq * ts * ts;

    //record the new setting 
    lp->current_rate = fosc;

    //integrated triangle LFO 
    lp->k_ = 2.0f * ts * frq;
    lp->nk_ = -2.0f * ts * frq;
    lp->psign_ = t;
    lp->nsign_ = -t;
    lp->sign_ = t;

    //Triangle LFO
    lp->ktri = frq / fs;

    //Sine LFO
    lp->ksin = PI * frq / fs;

    //Relaxation oscillator parameters
    float k = expf(-2.0f * fosc / fs);

    lp->rlx_k = k;
    lp->rlx_ik = 1.0f - k;

    //Exponential oscillator parameters
    k = expf(-2.0f * 1.3133f * fosc / fs);
    lp->exp_ik = k;
    lp->exp_k = 1.0f / k;

    if (lp->exp_x >= 1.0f) {
        lp->exp_x = lp->exp_k;
    } else {
        lp->exp_x = lp->exp_ik;
    }

    lp->exp_min = 1.0f / E;
    lp->exp_max = 1.0f + 1.0f / E;

    if (lp->exp_sv < lp->exp_min) {
        lp->exp_sv = lp->exp_min;
    }
    if (lp->exp_sv > lp->exp_max) {
        lp->exp_sv = lp->exp_max;
    }

}

//
// Integrated Triangle LFO shape 
//   The LFO shape appears sinusoidal but the derivative is 
//   a triangle wave.
//   This produces a linear (triangle) pitch shift when applied to 
//   a modulated delay line.
//

float run_integrated_triangle_lfo(lfoparams *lp) {
    //Wait for startup delay to begin at requested phase
    if (lp->startup_delay > 0) {
        lp->startup_delay -= 1;
        lp->lfo = 0.0f;

        return 0.0f;
    }

    lp->x += lp->sign;

    if (lp->x >= lp->k) {
        lp->sign = lp->nsign_;

        //Reset oscillator if rate change
        lp->x = lp->k_;
        lp->k = lp->k_;
        lp->nk = lp->nk_;

    } else if (lp->x <= lp->nk) {
        lp->sign = lp->psign_;

        //Reset oscillator if rate change
        lp->x = lp->nk_;
        lp->k = lp->k_;
        lp->nk = lp->nk_;
    }

    lp->lfo += lp->x;

    if (lp->lfo > 1.0f) {
        //fprintf(stderr,"%f\n",sum);
        lp->lfo = 1.0f;
    }
    if (lp->lfo < 0.0f) {
        //fprintf(stderr,"%f\n",sum);
        lp->lfo = 0.0f;
    }

    return lp->lfo;
}

// Triangle wave oscillator
float run_triangle_lfo(lfoparams *lp) {
    lp->trilfo += lp->ktri * lp->trisign;

    if (lp->trilfo >= 1.0f) {
        lp->trisign = -1.0f;
    }
    if (lp->trilfo <= 0.0f) {
        lp->trisign = 1.0f;
    }

    return lp->trilfo;
}

//Sine oscillator
float run_sine_lfo(lfoparams *lp) {
    lp->sin_part += lp->cos_part * lp->ksin;
    lp->cos_part -= lp->sin_part * lp->ksin;

    return 0.5f * (1.0f + lp->cos_part);
}


//
// RC Relaxation oscillator
// Simple 1rst-order lowpass filter driven by rlx_sign variable
// rlx_sign is toggled between 1.0 + 1/e and -1/e and
// toggles when LPF output reaches thresholds 0.0 and 1.0
// This relationship makes frequency computation simple
// as the oscillation period is simply 2*Tau (2RC)
//
// The waveshape is the typical RC filter (1-e^-t/Tau) charge
// shape followed by the RC filter exponential discharge shape.
// It appears roughly triangular.
// 
// This could be abused by adjusting rlx_min and rlx_max values
// to asymmetrical values, or closer to 1.0 to flatten the top more.
// 
float run_rlx_lfo(lfoparams *lp) {
    //Simple 1rst-order lowpass filter
    lp->rlx_lfo = lp->rlx_sign * lp->rlx_ik + lp->rlx_k * lp->rlx_lfo;

    //Schmitt trigger logic
    if (lp->rlx_lfo >= 1.0f) {
        lp->rlx_sign = lp->rlx_min;
    } else if (lp->rlx_lfo <= 0.0f) {
        lp->rlx_sign = lp->rlx_max;
    }

    //Output
    return lp->rlx_lfo;
}

float run_exp_lfo(lfoparams *lp) {
    //Exponential oscillator parameters
    lp->exp_sv *= lp->exp_x;

    if (lp->exp_sv >= lp->exp_max) {
        lp->exp_x = lp->exp_ik;
    } else if (lp->exp_sv <= lp->exp_min) {
        lp->exp_x = lp->exp_k;
    }

    return lp->exp_sv - lp->exp_min;
}

float run_lfo(lfoparams *lp) {
    float lfo_out;

    switch (lp->lfo_type) {
        case INT_TRI: //integrated triangle
            lfo_out = run_integrated_triangle_lfo(lp);
            break;

        case TRI: //triangle 
            lfo_out = run_triangle_lfo(lp);
            break;

        case SINE: //sine
            lfo_out = run_sine_lfo(lp);
            break;

        case SQUARE: //click-less square (compressed sine)
            lfo_out = run_sine_lfo(lp) - 0.5f;

            //Amplify and soft-clip the sine wave
            if (lfo_out > 0.0f) {
                lfo_out *= 1.0f / (1.0f + 30.0f * lfo_out);
            } else {
                lfo_out *= 1.0f / (1.0f - 30.0f * lfo_out);
            }

            lfo_out *= 16.0f;
            lfo_out += 0.5f;

            break;

        case EXP: //exponential
            lfo_out = run_exp_lfo(lp);
            break;

        case RELAX: //RC relaxation oscillator
            lfo_out = run_rlx_lfo(lp);
            break;

        case HYPER: //smooth bottom, triangular top
            lfo_out = run_integrated_triangle_lfo(lp);
            lfo_out = 1.0f - fabs(lfo_out - 0.5f);
            break;

        case HYPER_SINE:  //Sine bottom, triangular top
            lfo_out = run_sine_lfo(lp);
            lfo_out = 1.0f - fabs(lfo_out - 0.5f);
            break;

        default:
            lfo_out = run_integrated_triangle_lfo(lp);
            break;
    }

    return lfo_out;
}

void get_lfo_name(lfo_mode type, char *outstring) {
    // TODO: WTF? Return a string constant
    for (int i = 0; i < 30; i++) {
        outstring[i] = '\0';
    }

    switch (type) {
        case INT_TRI: //integrated triangle
            sprintf(outstring, "INTEGRATED TRIANGLE");
            break;

        case TRI: //triangle 
            sprintf(outstring, "TRIANGLE");
            break;

        case SINE: //sine
            sprintf(outstring, "SINE");
            break;

        case SQUARE: //click-less square (compressed sine)
            sprintf(outstring, "SQUARE");
            break;

        case EXP: //exponential
            sprintf(outstring, "EXPONENTIAL");
            break;

        case RELAX: //RC relaxation oscillator
            sprintf(outstring, "RC RELAXATION");
            break;

        case HYPER: //RC relaxation oscillator
            sprintf(outstring, "HYPER");
            break;

        case HYPER_SINE: //RC relaxation oscillator
            sprintf(outstring, "HYPER_SINE");
            break;

        default:
            sprintf(outstring, "DEFAULT: INTEGRATED TRIANGLE");
            break;
    }


}

void set_lfo_type(lfoparams *lp, lfo_mode type) {
    lp->lfo_type = type;
}
