/*
*   ADSR
*/
#include    <stdlib.h>

#include <cmath>

#include    "adsr.h"

adsr *make_adsr(float fs, int N) {
    adsr *ad;

    ad = (adsr *) malloc(sizeof(adsr));
    // System properties
    ad->fs = fs;
    ad->ifs = 1.0f / ad->fs;          //Sampling frequency and its inverse
    ad->NS = N;                   //Number of samples in loop

    // Parameters
    adsr_set_amplitude(ad, 1.0f);
    adsr_set_velocity(ad, 1.1414f);

    adsr_set_attack(ad, 50.0f);    //A
    adsr_set_decay(ad, 25.0f);     //D
    adsr_set_sustain(ad, 0.5f);    //S
    adsr_set_release(ad, 100.0f);  //R
    adsr_set_trigger_state(ad, false);

    // State variables
    adsr_set_trigger_timeout(ad, 0.25f);     //default to 250ms hold
    ad->state = ADSR_STATE_RELEASE;      //ADSR State machine current state
    ad->sv = 0.0f;

    return ad;
}

void adsr_set_amplitude(adsr *ad, float a) {
    ad->amplitude = a;
    ad->pk = ad->velocity * a;
    ad->sus = ad->amplitude * ad->sustain;
}

void adsr_set_velocity(adsr *ad, float v) {
    ad->velocity = v;
    ad->pk = ad->amplitude * v;
}

void adsr_set_attack(adsr *ad, float a) {
    float ms = a * 0.001f;

    //Limit 0.1ms to 1 second
    if (ms < 0.0001f) {
        ms = 0.0001f;
    } else if (ms > 1.0f) {
        ms = 1.0f;
    }

    ad->atk = 1.0f - expf(-ad->ifs / ms);
}

void adsr_set_decay(adsr *ad, float d) {
    float ms = d * 0.001f;

    //Limit 1ms to 2 seconds
    if (ms < 0.001) {
        ms = 0.001;
    } else if (ms > 2.0) {
        ms = 2.0;
    }

    ad->dcy = 1.0f - expf(-ad->ifs / ms);
}

void adsr_set_sustain(adsr *ad, float s) {
    float sus = s;

    //limit between 0 to 1
    if (sus > 1.0f) {
        sus = 1.0f;
    } else if (sus < 0.0f) {
        sus = 0.0f;
    }

    ad->sustain = sus;
    ad->sus = ad->amplitude * sus;
}

void adsr_set_release(adsr *ad, float r) {
    float ms = r * 0.001f;

    //Limit 1ms to 2 seconds
    if (ms < 0.001f) {
        ms = 0.001f;
    } else if (ms > 2.0f) {
        ms = 2.0f;
    }

    ad->rls = expf(-ad->ifs / ms);
}

//Set or clear trigger
void adsr_set_trigger_state(adsr *ad, bool t) {
    if (t) {
        if (ad->trig_timeout > 0) {
            ad->trig_timer = ad->trig_timeout;
        } else {
            ad->trig_timer = 2;
        }
        ad->state = ADSR_STATE_ATTACK;
    } else {
        ad->trig_timer = 0;
        ad->state = ADSR_STATE_RELEASE;
    }
}

void adsr_set_trigger_timeout(adsr *ad, float t) {
    int to = (int) (t * ad->fs);

    ad->trig_timeout = to;
}

void adsr_tick_n(adsr *ad, float *output) {
    for (int i = 0; i < ad->NS; i++) {
        if (ad->trig_timer < 1) {
            ad->state = ADSR_STATE_RELEASE;
        } else if (ad->trig_timeout > 0) {
            ad->trig_timer--;
        }

        //if trig_timeout < 0 then sustain forever
        if ((ad->state == ADSR_STATE_RELEASE) && (ad->trig_timer > 0)) {
            ad->state = ADSR_STATE_ATTACK;
        }

        switch (ad->state) {
            case ADSR_STATE_RELEASE:
                ad->sv *= ad->rls;
                break;

            case ADSR_STATE_ATTACK:
                if (ad->sv >= ad->amplitude) {
                    ad->sv = ad->amplitude;
                    ad->state = ADSR_STATE_DECAY;
                } else {
                    ad->sv = ad->sv + ad->atk * (ad->pk - ad->sv);
                }
                break;

            case ADSR_STATE_DECAY:
                ad->sv = ad->sv + ad->dcy * (ad->sus - ad->sv);
                if (ad->trig_timer == 0) {
                    ad->state = ADSR_STATE_RELEASE;
                }  //Redundant but paranoia safety
                break;

            default:
                ad->state = ADSR_STATE_RELEASE;  //more paranoia
                break;
        }
        output[i] = ad->sv;
    }
}
