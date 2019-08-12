//Band-limited interpolated delay line.  Uses 8th order butterworth IIR filter
//emulating combination of anti-aliasing & reconstruction filters from an 
//analog BBD typical circuit.

#include <stdlib.h>

#include <cmath>

#include "flange.h"

void tflanger_resize_delay(tflanger *cthis, float maxTime) {
    //printf("Updating delay line size.\n");
    size_t nSize = lrint(maxTime * cthis->fS + 1.0f);
    float *temp = (float *) malloc(nSize * sizeof(float));

    if (temp == nullptr) {
        // TODO: error handling
        // fprintf(stderr, "Cannot allocate any more memory.  Requested increase of delay time declined");
        return;
    }

    for (size_t i = 0; i < nSize; i++) {
        temp[i] = 0.0;
    }

    for (size_t i = 0; i < cthis->maxDly; i++) {
        temp[i] = cthis->dlyLine[i];
    }

    cthis->maxDly = nSize - 1;
    cthis->maxT = maxTime;

    float *dead = cthis->dlyLine;
    cthis->dlyLine = temp;

    free(dead);
}

tflanger *tflanger_init(float maxTime, float fSampleRate) {
    tflanger *cthis;

    cthis = (tflanger *) malloc(sizeof(tflanger));

    cthis->maxT = maxTime;
    cthis->fS = fSampleRate;

    const size_t nSize = lrint(cthis->maxT * fSampleRate + 1.0f);

    cthis->dlyLine = (float *) malloc(nSize * sizeof(float));

    for (size_t i = 0; i < nSize; i++) {
        cthis->dlyLine[i] = 0.0f;
    }

    cthis->maxDly = nSize - 1;
    //printf("maxT = %f\tmaxDly = %d\n", maxT, maxDly);

    //Delay line pointer
    cthis->dlyWrite = 0;

    //LFO
    cthis->lfopar = init_lfo(1.0f, cthis->fS, 0.0f);
    cthis->maxLfoRate = 100.0f;

    //Envelope detector
    cthis->envelope = (fparams *) malloc(sizeof(fparams));

    const float tc10 = 0.015f;  //About 10 Hz cutoff

    cthis->envelope->alpha = 1.0f / (tc10 * cthis->fS + 1.0f);
    cthis->envelope->ialpha = 1.0f - cthis->envelope->alpha;
    cthis->envelope->x1 = 0.0f;
    cthis->envelope->y1 = 0.0f;

    //Attack/Release settings
    cthis->attrel = (arparams *) malloc(sizeof(arparams));

    const float atkt = 0.05f;
    const float rlst = 0.5f;

    cthis->attrel->y1 = 0.0f;
    cthis->attrel->atk = 1.0f / (atkt * cthis->fS + 1.0f);
    cthis->attrel->iatk = 1.0f - cthis->attrel->atk;
    cthis->attrel->rls = 1.0f / (rlst * cthis->fS + 1.0f);
    cthis->attrel->irls = 1.0f - cthis->attrel->rls;

    //Bypass fader
    cthis->fader = (fparams *) malloc(sizeof(fparams));

    const float tcf = 0.25f;  //Fader time constant

    cthis->fader->alpha = 1.0f / (tcf * cthis->fS + 1.0f);
    cthis->fader->ialpha = 1.0f - cthis->fader->alpha;
    cthis->fader->x1 = 0.0f;
    cthis->fader->y1 = 0.0f;

    //Set some sane initial values
    const float ms = 0.001;

    tflanger_setLfoDepth(cthis, (9.0f * ms)); //lfo offset, input in seconds
    tflanger_setLfoWidth(cthis, (0.8f * ms)); //lfo deviation, input in seconds will
    //express peak-to-trough deviation
    tflanger_setLfoRate(cthis, 3.0f);    //Input in cycles per second
    tflanger_setLfoPhase(cthis, 0.0f);    //Phase offset, 0 to 2*pi, useful if
    //using 2 of these objects in stereo
    tflanger_setWetDry(cthis, 0.71f);     //give ratio wet, 0...1.0.
    //Dry = 1.0 - fracWet
    tflanger_setFeedBack(cthis, 0.0f);   //Regen

    //Bypass processing
    tflanger_setFinalGain(cthis, 1.0f);  // Set to 0.0 for bypass, set to >0.99 for norm
    cthis->trails = 1;                  //bypass with trails default

    //Envelope Detector Settings
    tflanger_setEnvelopeRateSkew(cthis, 0.0f);
    tflanger_setEnvelopeDepthSkew(cthis, 0.0f);
    tflanger_setEnvelopeWidthSkew(cthis, 0.0f);
    tflanger_setEnvelopeFbSkew(cthis, 0.0f);
    tflanger_setEnvelopeMixSkew(cthis, 0.0f);
    tflanger_setEnvelopeSensitivity(cthis, 1.0f);

    //Anti-aliasing filter
    const float fs = fSampleRate;
    const float f0 = 7200.0f;
    float *q;

    q = make_butterworth_coeffs(8);

    cthis->f[0] = make_biquad(LPF, fs, f0, q[0]);
    cthis->f[1] = make_biquad(LPF, fs, f0, q[1]);
    cthis->f[2] = make_biquad(LPF, fs, f0, q[2]);
    cthis->f[3] = make_biquad(LPF, fs, f0, q[3]);

    //Feedback tap to be run back through anti-alias filter with input
    cthis->regen = 0.0f;

    //Final computation
    tflanger_updateParams(cthis);

    return cthis;
}


void tflanger_destroy(tflanger *cthis) {
    free(cthis->dlyLine);
    free(cthis);
}

void get_wet_dry_mix(float fracWet, float *wet_out, float *dry_out) {
    float wet = 1.0f;
    float dry = 1.0f;
    float sign = 1.0f;

    if (fracWet >= 1.0f) {
        wet = 1.0f;
    } else if (fracWet < (-1.0f)) {
        wet = -1.0f; //allow add or subtract
    } else {
        wet = fracWet;
    }

    // TODO: use copysign
    if (wet < 0.0f) {
        sign = -1.0f;
    }

    float x = wet * wet;

    // make dry mix function
    float k = x * x;

    k *= k;
    dry = 1.0f - k;

    // make wet mix function
    k = 1 - x;
    k *= k;
    k *= k;
    wet = (1.0f - k) * sign;

    //keep dry at 0 phase
    *dry_out = dry;
    *wet_out = wet;
}

float potfunc1(float xn) {
    float x = 0.0f;
    int sign = 0;

    if (xn < 0.0f) {
        xn = -xn;
        sign = 1;
    }
    if (xn > 1.0f) {
        xn = 1.0f;
    }

    x = 2.0f * xn;

    if (x <= 1.0f) {
        xn = 0.5f * x * x;
    } else {
        xn = 1.0f - (1.0f - xn) * (2.0f - x);
    }

    // TODO: use copysign
    if (sign == 1) {
        xn = -xn;
    }

    return xn;
}


void tflanger_tick(tflanger *cthis, int nframes, float *samples, float *envelope) {
    if ((cthis->dry0 > 0.999f) && (cthis->trails == 0) && (cthis->outGain < 0.1f)) {
        return;
    }

    float dly;
    float dlyFloor;
    float fd;
    float ifd;

    float in; //current sample
    float envdet;

    long d = 0; //integer delay
    long d1 = 0; //integer delay - 1 

    for (int i = 0; i < nframes; i++) {
        //Envelope detector
        envdet = envelope[i];

        if (envdet > 1.0f) {
            envdet = 1.0f;
        } else if (envdet < 0.0f) {
            envdet = 0.0f;
        }

        envdet = tflanger_atkrls(cthis->attrel, envdet);

        // Bypass logic:  LPF set to slow time constant to fade input on and off
        // fades dry mix to unity in small steps at the same time.
        in = samples[i];
        tflanger_lpfilter(cthis->fader, cthis->outGain, 0);
        in *= cthis->fader->y1;

        if (cthis->outGain < 0.99f) {
            if (cthis->dry0 < 1.0f) {
                cthis->dry0 += 0.0001f;
            }  //-80dB steps to unity, about 1/4 second at 44.1 kHz
        }

        // Insert feedback
        in += cthis->regen;

        //anti-alias filter 
        for (int j = 0; j < 4; j++) {
            in = run_filter(in, cthis->f[j]);
        }

        //do ring buffer indexing
        if (++(cthis->dlyWrite) >= cthis->maxDly) {
            cthis->dlyWrite = 0;
        }

        cthis->dlyLine[cthis->dlyWrite] = in;


        //Apply to rate parameter
        if (cthis->rateskew != 0.0f) {
            float new_lfo_rate = cthis->rateskew * envdet + cthis->lfoRate;

            if (new_lfo_rate < 0.0f) {
                new_lfo_rate = 0.1f;
            }

            update_lfo(cthis->lfopar, new_lfo_rate, cthis->fS);
        } else if (cthis->lfoRate != cthis->lfopar->current_rate) {
            update_lfo(cthis->lfopar, cthis->lfoRate, cthis->fS);
        }

        //Apply to mix

        if (cthis->mixskew != 0.0f) {
            if (cthis->wet > 0.0f) {
                cthis->wet0 = cthis->mixskew * envdet + cthis->wet;
            } else {
                cthis->wet0 = cthis->wet - cthis->mixskew * envdet;
            }

            if (cthis->outGain >= 0.99f) {
                cthis->dry0 = cthis->dry - cthis->mixskew * envdet;
            }

            if (cthis->wet0 > 1.0f) {
                cthis->wet0 = 1.0f;
            }
            if (cthis->wet0 < -1.0f) {
                cthis->wet0 = -1.0f;
            }
            if (cthis->dry0 > 1.0f) {
                cthis->dry0 = 1.0f;
            }
            if (cthis->dry0 < -1.0f) {
                cthis->dry0 = -1.0f;
            }
        } else {
            cthis->wet0 = cthis->wet;

            if (cthis->outGain >= 0.99f) {
                cthis->dry0 = cthis->dry;
            }
        }

        //Apply to depth
        if (cthis->depthskew != 0.0f) {
            cthis->lfoDepth0 = cthis->depthskew * envdet + cthis->lfoDepth;

            if (cthis->lfoDepth0 < 0.0f) {
                cthis->lfoDepth0 = 0.0f;
            }
            if (cthis->lfoDepth0 >= cthis->maxT) {
                cthis->lfoDepth0 = cthis->maxT;
            }
        } else {
            cthis->lfoDepth0 = cthis->lfoDepth;
        }

        //Apply to width
        if (cthis->widthskew != 0.0f) {
            cthis->lfoWidth0 = cthis->widthskew * envdet + cthis->lfoWidth;

            if (cthis->lfoWidth0 < 0.0f) {
                cthis->lfoWidth0 = 0.0f;
            }

            if ((cthis->lfoWidth0 + cthis->lfoDepth0) >= cthis->maxT) {
                cthis->lfoWidth0 = cthis->maxT - cthis->lfoDepth0;
            }
        } else {
            cthis->lfoWidth0 = cthis->lfoWidth;
        }

        //Apply to feedback
        if (cthis->fbskew != 0.0f) {
            cthis->fb = cthis->fbskew * envdet + cthis->feedBack;

            if (cthis->fb >= 1.0f) {
                cthis->fb = 0.99f;    //less extreme limits
            }
            if (cthis->fb <= -1.0f) {
                cthis->fb = -0.99f;
            }
        } else {
            cthis->fb = cthis->feedBack;
        }

        //run LFO
        float lfo = run_lfo(cthis->lfopar);

        //calculate delay variables
        lfo *= cthis->lfoWidth0;
        dly = cthis->fS * (cthis->lfoDepth0 + lfo);
        dlyFloor = floor(dly);
        fd = (dly - dlyFloor);
        ifd = 1.0f - fd;

        //set up high delay tap
        d = lrint(dlyFloor);
        d1 = d + 1;
        d = cthis->dlyWrite - d;
        d1 = cthis->dlyWrite - d1;

        if (d < 0) {
            d += cthis->maxDly;
        }
        if (d1 < 0) {
            d1 += cthis->maxDly;
        }

        //output signal
        cthis->regen = fd * cthis->dlyLine[d1] + ifd * cthis->dlyLine[d];

        if (cthis->trails == 0) {
            samples[i] = cthis->dry0 * samples[i] + cthis->wet0 * cthis->fader->y1 * cthis->regen;
        } else {
            samples[i] = cthis->dry0 * samples[i] + cthis->wet0 * cthis->regen;
        }

        cthis->regen *= cthis->fb;
    } //for(i...
}


float tflanger_lpfilter(fparams *cthis, float data, char mode) {
    const float y0 = cthis->alpha * data + cthis->ialpha * cthis->y1;

    cthis->y1 = y0;

    // TODO: use enum/constant for mode
    if (mode == 0) {
        return y0;
    } else {
        return (0.5f * (y0 + data));   //low-pass shelf where high freqs shelve to 1/2 value of low freqs
    }
}

float tflanger_atkrls(arparams *cthis, float data) {
    float y0;

    if (data > cthis->y1) {
        y0 = cthis->atk * data + cthis->iatk * cthis->y1;  //Attack
    } else {
        y0 = cthis->rls * data + cthis->irls * cthis->y1;  //Release
    }

    cthis->y1 = y0;

    return y0;
}

void tflanger_updateParams(tflanger *cthis) {
    if ((cthis->lfoWidth + cthis->lfoDepth) >= cthis->maxT) {
        cthis->lfoWidth = cthis->maxT - cthis->lfoDepth;
    }

    /*
    cthis->lfoWidth0 = cthis->lfoWidth;
    cthis->lfoDepth0 = cthis->lfoDepth;
    cthis->wet0 = cthis->wet;
    cthis->dry0 = cthis->dry;
    cthis->fb = cthis->feedBack;
    */

    //printf("\n\nLFO Depth = %lf\nLFO Width = %lf\nLFO Rate = %lf\n%%Wet:  %lf\n%%Feedback: %lf\nOutput Level:  %lf\n", cthis->lfoDepth, cthis->lfoWidth, cthis->lfoRate, cthis->wet, cthis->feedBack, cthis->outGain);
}

void tflanger_setLfoDepth(tflanger *cthis, float lfoDepth_) {
    cthis->lfoDepth = fabs(lfoDepth_);

    //printf("d_ = %lf\tdepth = %lf\n\n", lfoDepth_, lfoDepth);

    if ((2.0f * cthis->lfoDepth) >= cthis->maxT) {
        tflanger_resize_delay(cthis, 2.0f * cthis->lfoDepth);
    }

    tflanger_updateParams(cthis);
    // printf("d_ = %lf\tdepth = %lf\n\n", lfoDepth_, lfoDepth);
}

void tflanger_setLfoWidth(tflanger *cthis, float lfoWidth_) {
    cthis->lfoWidth = fabs(lfoWidth_);

    if ((cthis->lfoWidth + cthis->lfoDepth) >= cthis->maxT) {
        tflanger_resize_delay(cthis, cthis->lfoWidth + cthis->lfoDepth);
    }

    tflanger_updateParams(cthis);
}

void tflanger_setLfoRate(tflanger *cthis, float lfoRate_) {
    cthis->lfoRate = fabs(lfoRate_);

    if (cthis->lfoRate >= cthis->maxLfoRate) {
        cthis->lfoRate = cthis->maxLfoRate;
    }

    update_lfo(cthis->lfopar, cthis->lfoRate, cthis->fS);
}

void tflanger_setLfoPhase(tflanger *cthis, float lfoPhase_) {
    cthis->lfoPhase = lfoPhase_;
    cthis->lfoPhase = fmod(cthis->lfoPhase, 2.0f * PI);

    tflanger_updateParams(cthis);
}

void tflanger_setWetDry(tflanger *cthis, float fracWet) {
    float wet = 1.0f;
    float dry = 1.0f;

    get_wet_dry_mix(fracWet, &wet, &dry);

    //keep dry at 0 phase
    cthis->dry = dry;
    cthis->wet = wet;

    tflanger_updateParams(cthis);
}

void tflanger_setFeedBack(tflanger *cthis, float feedBack_) {
    cthis->feedBack = feedBack_;//potfunc1(feedBack_);//*(2.0 - fabs(feedBack_));

    if (cthis->feedBack >= 0.99f) {
        cthis->feedBack = 0.99f;
    } else if (cthis->feedBack <= -0.99f) {
        cthis->feedBack = -0.99f;
    }
}

void tflanger_setDamping(tflanger *cthis, float fdamp_) {
    float fdamp = fdamp_;

    //Q of the final stage comes out ot about 0.51 
    //so it is a relatively gentle roll-off
    float Q = cthis->f[3]->Q;

    if (fdamp < 40.0f) {
        fdamp = 40.0f;
    } else if (fdamp > 7200.0f) {
        fdamp = 7200.0f;
    }

    // 
    biquad_update_coeffs(LPF, cthis->f[3], cthis->fS, fdamp, Q);
}


void tflanger_setFinalGain(tflanger *cthis, float outGain_) {
    cthis->outGain = outGain_;

    //limit gain here to something sane
    if (cthis->outGain >= 1.42f) {
        cthis->outGain = 1.42f;  //1.42 is about +3dB
    }
    if (cthis->outGain < 0.0f) {
        cthis->outGain = 0.0f;
    }
}

void tflanger_setTrails(tflanger *cthis, char trails) {
    cthis->trails = trails;
}

void tflanger_setEnvelopeSensitivity(tflanger *cthis, float sns) {
    if (sns < 0.0f) {
        cthis->envelope_sensitivity = 0.0f;
    } else if (sns > 36.0f) {
        sns = 36.0f;
    } else { // TODO: sns never used. remove else?
        cthis->envelope_sensitivity = sns;
    }
}

void tflanger_setEnvelopeAttack(tflanger *cthis, float atk) {
    cthis->attrel->atk = 1.0f / (atk * cthis->fS + 1.0f);
    cthis->attrel->iatk = 1.0f - cthis->attrel->atk;
}

void tflanger_setEnvelopeRelease(tflanger *cthis, float rls) {
    cthis->attrel->rls = 1.0f / (rls * cthis->fS + 1.0f);
    cthis->attrel->irls = 1.0f - cthis->attrel->rls;
}

void tflanger_setEnvelopeRateSkew(tflanger *cthis, float skew) {
    cthis->rateskew = skew;

    if (skew > 100.0f) {
        cthis->rateskew = 100.0f;
    }
    if (skew < -100.0f) {
        cthis->rateskew = -100.0f;
    }
}

void tflanger_setEnvelopeDepthSkew(tflanger *cthis, float skew) {
    cthis->depthskew = skew;

    if (skew > 1.0f) {
        cthis->depthskew = 1.0f;
    }
    if (skew < -1.0f) {
        cthis->depthskew = -1.0f;
    }
}

void tflanger_setEnvelopeWidthSkew(tflanger *cthis, float skew) {
    cthis->widthskew = skew;

    if (skew > 1.0f) {
        cthis->widthskew = 1.0f;
    }
    if (skew < -1.0f) {
        cthis->widthskew = -1.0f;
    }
}

void tflanger_setEnvelopeFbSkew(tflanger *cthis, float skew) {
    cthis->fbskew = skew;

    if (skew > 1.0f) {
        cthis->fbskew = 1.0f;
    }
    if (skew < -1.0f) {
        cthis->fbskew = -1.0f;
    }
}


void tflanger_setEnvelopeMixSkew(tflanger *cthis, float skew) {
    cthis->mixskew = skew;

    if (skew > 1.0f) {
        cthis->mixskew = 1.0f;
    }
    if (skew < -1.0f) {
        cthis->mixskew = -1.0f;
    }
}

void tflanger_set_lfo_type(tflanger *cthis, unsigned int type) {
    set_lfo_type(cthis->lfopar, type);
}

void tflanger_setPreset(tflanger *cthis, unsigned int preset) {
    float ms = 0.001f;

    switch (preset) {
        // TODO: use constants/enum instead of hardcoded values
        case 0: //basic chorus
            //Set some sane initial values

            tflanger_setLfoDepth(cthis, (9.5f * ms)); //lfo offset, input in seconds
            tflanger_setLfoWidth(cthis, (1.2f * ms)); //lfo deviation, input in seconds will
            //express peak-to-trough deviation
            tflanger_setLfoRate(cthis, 2.6f);    //Input in cycles per second
            //LFO Rate synchronized to default delay preset time
            tflanger_setLfoPhase(cthis, 0.0f);    //Phase offset, 0 to 2*pi, useful if
            //using 2 of these objects in stereo
            tflanger_setWetDry(cthis, 0.5f);     //give ratio wet, 0...1.0.
            //Dry = 1.0 - fracWet
            tflanger_setFeedBack(cthis, 0.0f);   //Regen
            tflanger_setTrails(cthis, 0);
            break;

        case 1: //basic delay
            tflanger_setLfoDepth(cthis, 385.0f * ms);   //delay time
            tflanger_setLfoWidth(cthis, 0.25f * ms); //mild modulation
            tflanger_setFeedBack(cthis, 0.36f);
            tflanger_setWetDry(cthis, 0.45f);
            tflanger_setLfoRate(cthis, 1.3f);
            tflanger_setDamping(cthis, 2800.0f);  //Similar range to typical BBD delay pedal
            tflanger_setTrails(cthis, 1);
            break;

        case 2: //basic flanger 
            tflanger_setLfoDepth(cthis, 0.8f * ms);   //sweeps pretty high
            tflanger_setLfoWidth(cthis, 4.0f * ms); //mild modulation
            tflanger_setFeedBack(cthis, -0.3f);
            tflanger_setWetDry(cthis, -0.5f);
            tflanger_setLfoRate(cthis, 0.325f); //LFO Rate synchronized to default delay preset time
            tflanger_setTrails(cthis, 0);
            break;

        default:
            break;
    }
}
