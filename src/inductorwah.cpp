#include <cmath>

#include "inductorwah.h"

void commit_circuit_config(iwah_coeffs *cf, float fs) {
    //Useful constants
    const float RpRi = cf->Rp * cf->Ri / (cf->Rp + cf->Ri);
    const float f0 = 1.0f / (2.0f * PI * std::sqrt(cf->Lp * cf->Cf));
    const float w0 = 2.0f * PI * f0 / fs;
    const float Q = RpRi * std::sqrt(cf->Cf / cf->Lp);
    const float c = std::cos(w0);
    const float s = std::sin(w0);
    const float alpha = s / (2.0f * Q);

    //High Pass Biquad Coefficients
    cf->b0h = (1.0f + c) / 2.0f;
    cf->b1h = -(1.0f + c);
    cf->b2h = (1.0f + c) / 2.0f;

    //Band-pass biquad coefficients
    cf->b0b = Q * alpha;
    cf->b2b = -Q * alpha;
    cf->a0b = 1.0f + alpha;
    cf->a1b = -2.0f * c;
    cf->a2b = 1.0f - alpha;


    //1-pole high pass filter coefficients
    // H(z) = g * (1 - z^-1)/(1 - a1*z^-1)
    // Direct Form 1:
    //      y[n] = ghpf * ( x[n] - x[n-1] ) - a1p*y[n-1]

    cf->a1p = -expf(-1 / (cf->Ri * cf->Ci * fs));
    cf->ghpf = cf->gf * (1.0f - cf->a1p) * 0.5f;   //BJT forward gain worked in here to
    // save extra multiplications in
    // updating biquad coefficients

    //Distill all down to final biquad coefficients
    const float Gi = cf->Rs / (cf->Ri + cf->Rs);
    const float gbpf = 1.0f / (2.0f * PI * f0 * cf->Ri * cf->Cf);  //band-pass component equivalent gain

    //Final Biquad numerator coefficients
    cf->b0 = gbpf * cf->b0b + Gi * cf->a0b;
    cf->b1 = Gi * cf->a1b;
    cf->b2 = gbpf * cf->b2b + Gi * cf->a2b;

    //Constants to make denominator coefficients computation more efficient
    //in real-time
    cf->a0c = -cf->gf * cf->b0h;
    cf->a1c = -cf->gf * cf->b1h;
    cf->a2c = -cf->gf * cf->b2h;

    //denominator coefficients (initialization)
    cf->a0 = cf->a0b;
    cf->a1 = -cf->a1b; //negated because of addition used in
    cf->a2 = -cf->a2b; //biquad filter form
    cf->ax = 1.0f / cf->a0;

}

void iwah_circuit_preset(iwah_mode ckt, iwah_coeffs *cf, float fs) {
    //Useful multipliers for short-hand expression of component values
    const float k = 1000.0f;
    const float n = 1.0e-9f;

    //helper variables 
    float ro;
    float re;
    float Req;
    float ic;
    float beta;

    switch (ckt) {
        case GCB: //Dunlop GCB-95 variant 
            //Default values from Crybaby GCB-95 circuit
            cf->Rc = 22.0f * k; //BJT gain stage collector resistor
            cf->Rpot = 100.0f * k; //Pot resistance value
            cf->Rbias = 470.0f * k; //Typically 470k bias resistor shows up in parallel with output
            cf->Re = 390.0f; //BJT gain stage emitter resistor

            //Equivalent output resistance seen from BJT collector
            ro = cf->Rc * cf->Rpot / (cf->Rc + cf->Rpot);
            ro = ro * cf->Rbias / (ro + cf->Rbias);

            ic = 3.7f / cf->Rc;  //Typical bias current
            re = 0.025f / ic; //BJT gain stage equivalent internal emitter resistance
            // gm = Ic/Vt, re = 1/gm
            Req = cf->Re + re;

            cf->gf = -ro / Req;  //forward gain of transistor stage

            beta = 250.0f;  //BJT forward gain
            cf->re = beta * Req; //Resistance looking into BJT emitter

            cf->Cf = 10.0f * n;
            cf->Ci = 10.0f * n;
            cf->Rp = 33.0f * k * cf->re / (cf->re + 33.0f * k);
            cf->Lp = 0.5f;
            cf->Ri = 68.0f * k;
            cf->Rs = 1.5f * k;
            break;

        case VOX: //Vox V-847 voicing
            cf->Rc = 22.0f * k; //BJT gain stage collector resistor
            cf->Rpot = 100.0f * k; //Pot resistance value
            cf->Rbias = 470.0f * k; //Typically 470k bias resistor shows up in parallel with output
            cf->Re = 510.0f; //BJT gain stage emitter resistor

            //Equivalent output resistance seen from BJT collector
            ro = cf->Rc * cf->Rpot / (cf->Rc + cf->Rpot);
            ro = ro * cf->Rbias / (ro + cf->Rbias);

            ic = 3.7f / cf->Rc;  //Typical bias current
            re = 0.025f / ic; //BJT gain stage equivalent internal emitter resistance
            // gm = Ic/Vt, re = 1/gm
            Req = cf->Re + re;

            cf->gf = -ro / Req;  //forward gain of transistor stage

            beta = 650.0f;  //BJT forward gain
            cf->re = beta * Req; //Resistance looking into BJT emitter

            cf->Cf = 10.0f * n;
            cf->Ci = 10.0f * n;
            cf->Ri = 68.0f * k;
            cf->Rp = 33.0f * k;
            cf->Rp = cf->Rp * cf->Ri / (cf->Rp + cf->Ri);
            cf->Lp = 0.5f;
            cf->Rs = 1.5f * k;
            break;

        case DUNLOP: //Original Crybaby voicing
            cf->Rc = 22.0f * k; //BJT gain stage collector resistor
            cf->Rpot = 100.0f * k; //Pot resistance value
            cf->Rbias = 470.0f * k; //Typically 470k bias resistor shows up in parallel with output
            cf->Re = 470.0f; //BJT gain stage emitter resistor

            //Equivalent output resistance seen from BJT collector
            ro = cf->Rc * cf->Rpot / (cf->Rc + cf->Rpot);
            ro = ro * cf->Rbias / (ro + cf->Rbias);

            ic = 3.7f / cf->Rc;  //Typical bias current
            re = 0.025f / ic; //BJT gain stage equivalent internal emitter resistance
            // gm = Ic/Vt, re = 1/gm
            Req = cf->Re + re;

            cf->gf = -ro / Req;  //forward gain of transistor stage

            beta = 250.0f;  //BJT forward gain
            cf->re = beta * Req; //Resistance looking into BJT emitter

            cf->Cf = 10.0f * n;
            cf->Ci = 10.0f * n;
            cf->Rp = 33.0f * k * cf->re / (cf->re + 33.0f * k);
            cf->Lp = 0.66f;
            cf->Ri = 68.0f * k;
            cf->Rs = 1.5f * k;
            break;

        case MCCOY: //Clyde McCoy voicing
            cf->Rc = 22.0f * k; //BJT gain stage collector resistor
            cf->Rpot = 100.0f * k; //Pot resistance value
            cf->Rbias = 470.0f * k; //Typically 470k bias resistor shows up in parallel with output
            cf->Re = 470.0f; //BJT gain stage emitter resistor

            //Equivalent output resistance seen from BJT collector
            ro = cf->Rc * cf->Rpot / (cf->Rc + cf->Rpot);
            ro = ro * cf->Rbias / (ro + cf->Rbias);

            ic = 3.7f / cf->Rc;  //Typical bias current
            re = 0.025f / ic; //BJT gain stage equivalent internal emitter resistance
            // gm = Ic/Vt, re = 1/gm
            Req = cf->Re + re;

            cf->gf = -ro / Req;  //forward gain of transistor stage

            beta = 250.0f;  //BJT forward gain
            cf->re = beta * Req; //Resistance looking into BJT emitter

            cf->Cf = 10.0f * n;
            cf->Ci = 10.0f * n;
            cf->Rp = 100.0f * k * cf->re / (cf->re + 100.0f * k);
            cf->Lp = 0.5f;
            cf->Ri = 68.0f * k;
            cf->Rs = 1.5f * k;
            break;

        case VOCAL: //Vox with a little more vocal voicing that also cuts through distortion
            cf->Rc = 22.0f * k; //BJT gain stage collector resistor
            cf->Rpot = 100.0f * k; //Pot resistance value
            cf->Rbias = 470.0f * k; //Typically 470k bias resistor shows up in parallel with output
            cf->Re = 510.0f; //BJT gain stage emitter resistor

            //Equivalent output resistance seen from BJT collector
            ro = cf->Rc * cf->Rpot / (cf->Rc + cf->Rpot);
            ro = ro * cf->Rbias / (ro + cf->Rbias);

            ic = 3.7f / cf->Rc;  //Typical bias current
            re = 0.025f / ic; //BJT gain stage equivalent internal emitter resistance
            // gm = Ic/Vt, re = 1/gm
            Req = cf->Re + re;

            cf->gf = -ro / Req;  //forward gain of transistor stage

            beta = 250.0f;  //BJT forward gain
            cf->re = beta * Req; //Resistance looking into BJT emitter

            cf->Cf = 15.0f * n;
            cf->Ci = 8.0f * n;
            cf->Rp = 47.0f * k * cf->re / (cf->re + 47.0f * k);
            cf->Lp = 0.5f;
            cf->Ri = 68.0f * k;
            cf->Rs = 800.0f;
            break;

        case EXTREME: //Vox with a little more vocal voicing that also cuts through distortion
            cf->Rc = 22.0f * k; //BJT gain stage collector resistor
            cf->Rpot = 100.0f * k; //Pot resistance value
            cf->Rbias = 470.0f * k; //Typically 470k bias resistor shows up in parallel with output
            cf->Re = 150.0f; //BJT gain stage emitter resistor

            //Equivalent output resistance seen from BJT collector
            ro = cf->Rc * cf->Rpot / (cf->Rc + cf->Rpot);
            ro = ro * cf->Rbias / (ro + cf->Rbias);

            ic = 3.7f / cf->Rc;  //Typical bias current
            re = 0.025f / ic; //BJT gain stage equivalent internal emitter resistance
            // gm = Ic/Vt, re = 1/gm
            Req = cf->Re + re;

            cf->gf = -ro / Req;  //forward gain of transistor stage

            beta = 1200.0f;  //BJT forward gain
            cf->re = beta * Req; //Resistance looking into BJT emitter

            cf->Cf = 10.0f * n;
            cf->Ci = 47.0f * n;
            cf->Rp = 150.0f * k;
            cf->Lp = 0.5;
            cf->Ri = 220.0f * k;
            cf->Rs = 100;
            break;

        default:
            break;
    }

    //Compute DSP filter coefficients
    commit_circuit_config(cf, fs);
}

//Initialize filter state variables
void zero_state_variables(iwah_coeffs *cf) {
    //biquad state variables
    cf->y1 = 0.0f;
    cf->y2 = 0.0f;
    cf->x1 = 0.0f;
    cf->x2 = 0.0f;

    //First order high-pass filter state variables
    cf->yh1 = 0.0f;
    cf->xh1 = 0.0f;
}

iwah_coeffs *make_iwah(float fs) {
    iwah_coeffs *cf;

    //First malloc the struct
    cf = (iwah_coeffs *) malloc(sizeof(iwah_coeffs));

    //Default to GCB-95 circuit
    iwah_circuit_preset(GCB, cf, fs);

    //Initialize everything to 0.0
    zero_state_variables(cf);

    //setup as bypassed by default 
    cf->bypass = true;

    return cf;
}

inline float sqr(float x) {
    return x * x;
}

inline float clip1(float x) {
    const float thrs = 0.8f;
    const float nthrs = -0.72f;
    const float f = 1.25f;

    //Hard limiting
    if (x >= 1.2f) {
        x = 1.2f;
    }
    if (x <= -1.12f) {
        x = -1.12f;
    }

    //Soft clipping
    if (x > thrs) {
        x -= f * sqr(x - thrs);
    }
    if (x < nthrs) {
        x += f * sqr(x - nthrs);
    }

    return x;
}

inline float iwah_tick(float x, float gp, iwah_coeffs *cf) {
    //variable gp is the pot gain, nominal range 0.0 to 1.0
    //although this can be abused for extended range.
    //A value less than zero would make the filter go unstable
    //and put a lot of NaNs in your audio output
    if (gp < 0.0f) {
        gp = 0.0f;
    }

    //The magic numbers below approximate frequency warping characteristic
    const float gw = 4.6f - 18.4f / (4.0f + gp);

    //Update Biquad coefficients
    cf->a0 = cf->a0b + gw * cf->a0c;
    cf->a1 = -(cf->a1b + gw * cf->a1c);
    cf->a2 = -(cf->a2b + gw * cf->a2c);
    cf->ax = 1.0f / cf->a0;

    //run it through the 1-pole HPF and gain first
    const float hpf = cf->ghpf * (x - cf->xh1) - cf->a1p * cf->yh1;

    cf->xh1 = x;
    cf->yh1 = hpf;

    //Apply modulated biquad
    float y0 = (cf->b0 * hpf + cf->b1 * cf->x1 + cf->b2 * cf->x2
               + (cf->a1 * cf->y1 + cf->a2 * cf->y2)) * cf->ax;

    y0 *= cf->ax;

    const float out = clip1(y0);

    y0 = 0.95f * y0 + 0.05f * out; //Let a little harmonic distortion feed back into the filter loop
    cf->x2 = cf->x1;
    cf->x1 = hpf;
    cf->y2 = cf->y1;
    cf->y1 = y0;

    //And ship it!
    return out;

}

void iwah_tick_n(iwah_coeffs *cf, float *x, float *gp, int n) {
    if (cf->bypass) {
        zero_state_variables(cf);
        return;
    }

    for (int i = 0; i < n; i++) {
        x[i] = iwah_tick(x[i], gp[i], cf);
    }
}

void iwah_bypass(iwah_coeffs *cf, bool state) {
    cf->bypass = state;
}
