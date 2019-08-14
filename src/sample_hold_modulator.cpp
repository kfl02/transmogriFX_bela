//
//   Sample/hold modulator program
//
// -----------------------------------------------------------------------------
// This is free and unencumbered software released into the public domain.
//
// For more information, please refer to <http://unlicense.org/>
// terms re-iterated in UNLICENSE.txt
// -----------------------------------------------------------------------------
//

#include    <stdlib.h>

#include    "sample_hold_modulator.h"

inline float sqr(float x) {
    return x * x;
}

inline float soft_clip(float xn_) {
    float xn = 2.0f * xn_ - 1.0f;

    if (xn > 1.0f) {
        xn = 1.0f;
    } else if (xn < -1.0f) {
        xn = -1.0f;
    } else if (xn < 0.0f) {
        xn = sqr(xn + 1.0f) - 1.0f;
    } else {
        xn = 1.0f - sqr(1.0f - xn);
    }

    return 0.5f * (xn + 1.0f);
}

sh_mod *make_sample_hold(float fs, int N) {
    sh_mod *sh;

    sh = (sh_mod *) malloc(sizeof(sh_mod));

    sh->ad = make_adsr(fs, N);
    sh->max_types = SH_MAX_TYPES;
    sh->adsr_env = (float *) malloc(sizeof(float) * N);
    sh->lfo = (float *) malloc(sizeof(float) * N);
    sh->hfo = (float *) malloc(sizeof(float) * N);
    sh->swave = (float *) malloc(sizeof(float) * N);
    sh->sequence = (float *) malloc(sizeof(float) * MAX_SEQ_STEPS);
    sh->sclk = (bool *) malloc(sizeof(bool) * N);

    for (int i = 0; i < N; i++) {
        sh->adsr_env[i] = 0.0f;
        sh->lfo[i] = 0.0f;
        sh->hfo[i] = 0.0f;
        sh->swave[i] = 0.0f;
        sh->sclk[i] = false;
    }

    sh->fs = fs;
    sh->ifs = 1.0f / sh->fs;
    sh->NS = N;

    const float rate = 6.0f;
    const float hfr = 4.7f * rate;
    const float Trate = 1000.0f / rate;

    //Initialize ramp rates
    adsr_set_attack(sh->ad, Trate / 8.0f);
    adsr_set_decay(sh->ad, Trate / 8.0f);
    adsr_set_sustain(sh->ad, 0.36f);
    adsr_set_release(sh->ad, Trate / 4.0f);
    adsr_set_trigger_timeout(sh->ad, 0.5f / rate);
    sh->en_adsr = true;

    //Initialize ramp rates
    sh->dt_lfo = sh->ifs * rate;
    sh->dt_hfo = sh->ifs * hfr;
    sh->steepness = 3.0f;  //how abruptly it transitions from one step to the next
    sh->mode = SH_RAND;
    sh->n_seq_steps = 6;

    for (int i = 0; i < sh->n_seq_steps; i++) {
        sh->sequence[i] = ((float) rand()) / ((float) RAND_MAX);
        //sh->sequence[i] = ((float) i)/((float) (sh->n_seq_steps -1) );
    }

    //Initialize state variables
    sh->seq_step = 0;
    sh->ramp_lfo = 0.0f;
    sh->ramp_hfo = 0.0f;
    sh->last_hfo_sample = 0.0f;
    sh->last_wave = 0.0f;
    sh->dl = 0.0f;

    return sh;
}

void run_modulator(sh_mod *sh) {
    for (int i = 0; i < sh->NS; i++) {
        sh->ramp_lfo += sh->dt_lfo;

        if (sh->mode == SH_RAND) {
            sh->ramp_hfo = ((float) rand()) / ((float) RAND_MAX);
        } else if (sh->mode == SH_RAMP) {
            sh->ramp_hfo += sh->dt_hfo;
        } else if (sh->mode == SH_SEQ) {
            sh->ramp_hfo = sh->sequence[sh->seq_step];
        }

        if (sh->ramp_lfo >= 1.0f) {
            sh->ramp_lfo = 0.0f;
            sh->sclk[i] = true;
            sh->seq_step += 1;  //go to next value sequence
            if (sh->seq_step >= sh->n_seq_steps) {
                sh->seq_step = 0;  //reset sequencer
            }
        } else {
            sh->sclk[i] = false;
        }

        if (sh->ramp_hfo > 1.0f) {
            sh->ramp_hfo = 1.0f;
        } else if (sh->ramp_hfo < 0.0f) {
            sh->ramp_hfo = 0.0f;
        }

        sh->lfo[i] = soft_clip(sh->steepness * sh->ramp_lfo);
        sh->hfo[i] = sh->ramp_hfo;
    }
}

void run_sample_hold(sh_mod *sh, float *output) {
    run_modulator(sh);
    adsr_tick_n(sh->ad, sh->adsr_env);

    for (int i = 0; i < sh->NS; i++) {
        if (sh->sclk[i]) {
            sh->dl = (sh->hfo[i] - sh->last_hfo_sample);
            sh->last_wave = sh->last_hfo_sample;
            sh->last_hfo_sample = sh->hfo[i];

            adsr_set_trigger_state(sh->ad, true);
            adsr_set_amplitude(sh->ad, sh->last_hfo_sample);
        }

        sh->swave[i] = sh->dl * sh->lfo[i] + sh->last_wave;

        if (sh->en_adsr) {
            output[i] = sh->adsr_env[i];
            //output[i] = sh->swave[i]*sh->adsr_env[i];
            //output[i] = sh->last_hfo_sample*sh->adsr_env[i];
            //output[i] = sh->last_hfo_sample;
        } else {
            output[i] = sh->swave[i];
        }
    }
}

bool sample_hold_set_active(sh_mod *sh, bool act) {
    //toggle state if act is true
    if (act) {
        sh->en_adsr = !sh->en_adsr;
    } else {
        //always deactivate if false
        sh->en_adsr = false;
    }

    //return the state evaluated
    return sh->en_adsr;
}

void sample_hold_set_rate(sh_mod *sh, float rate) {
    const float hfr = 4.7f * rate;
    const float Trate = 1000.0f / rate;

    //Initialize ramp rates
    adsr_set_attack(sh->ad, Trate / 10.0f);
    adsr_set_decay(sh->ad, Trate / 20.0f);
    adsr_set_sustain(sh->ad, 0.75f);
    adsr_set_release(sh->ad, Trate / 16.0f);
    adsr_set_trigger_timeout(sh->ad, 0.5f / rate);
    sh->en_adsr = true;

    //Initialize ramp rates
    sh->dt_lfo = sh->ifs * rate;
    sh->dt_hfo = sh->ifs * hfr;
}

void sample_hold_set_type(sh_mod *sh, sh_mod_mode type) {
    sh_mod_mode t = type;

    if (t >= SH_MAX_TYPES) {
        t = SH_SEQ;
    } else if (t < 0) {
        t = SH_RAND;
    }
     sh->mode = t;
}
