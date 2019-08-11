//
// KOT style tonestack
//
#include "kot_tonestack.h"

void kotstack_compute_filter_coeffs(kot_stack *ks) {
    // Local variables to make following formulas easier
    // Tone pot ratios
    const float apos = ks->tone_pot_pos;
    const float bpos = 1.0f - apos;  // Alternate pot leg resistance ratio
    //printf("apos= %f, bpos= %f\n", apos, bpos);

    // Boost pot ratio
    const float bstpos = ks->boost_pot_pos;

    // Sampling rate
    const float fs = ks->fs;
    const float kk = 2.0f * fs;
    //printf("fs= %f, kk=%f\n", fs, kk);

    // Circuit params
    const float rtone = ks->rtone;
    const float rboost = ks->rboost;
    const float ri = ks->ri;
    const float rc = ks->rc;
    const float ro = ks->ro;
    const float c1 = ks->c1;
    const float c2 = ks->c2;

    // Components as labeled in derivation
    const float ra = ri + apos * rtone;
    const float rb = bpos * rtone + rtone / fs; // This sets one of the zeros, so it needs to avoid going to zero
    const float rd = rboost * bstpos;
    //printf("ra= %f, rb= %f, rd= %f\n", ra, rb, rd);

    // Substitutions as labled in derivation
    const float a0 = ro + rc;
    const float a1 = (ro + rd) * rc * c2 + ro * rd * c2;
    const float b1 = (ro + rd) * c2;
    //printf("a0= %f, a1= %f, b1= %f\n", a0, a1, b1);

    const float N0 = a0;
    const float N1 = rb * c1 * a0 + a1;
    const float N2 = rb * c1 * a1;
    const float D0 = 1.0;
    const float D1 = rb * c1 + c1 * a0 + b1;
    const float D2 = rb * c1 * b1 + c1 * a1;
    //printf("N0= %f, N1= %f, N2= %f, D0= %f, D1= %f, D2= %f\n", N0, N1, N2, D0, D1, D2);

    // Analog Filter consolidated coefficients

    // Stage 1 filter (order 2)
    // vxvi = gvx*(s*s + A1*s + A0)/(s*s + B1*s + B0)
    const float gvx = N2 / (N2 + ra * D2);
    const float A0 = N0 / N2;
    const float A1 = N1 / N2;
    const float B0 = (N0 + ra * D0) / (N2 + ra * D2);
    const float B1 = (N1 + ra * D1) / (N2 + ra * D2);
    //printf("A0= %f, A1= %f, B0= %f, B1= %f, gvx= %f\n", A0, A1, B0, B1, gvx);

    // Stage 2 filter (order 1)
    // vovx = gvo*(s + X0)/(s + Y0);
    const float gvo = (c2 * ro * rd) / (c2 * (ro * rd + (ro + rd) * rc));
    const float X0 = 1 / (c2 * rd);
    const float Y0 = (rc + ro) / (c2 * (ro * rd + (ro + rd) * rc));

    // Bilinear transform computations

    // Stage 1 filter (order 2)
    // vxviblt = gdsp*(A2blt*z2 + A1blt*z1 + 1) / ( B2blt*z2 + B1blt*z1 + 1)
    const float gblt = ((A0 + kk * A1 + kk * kk) / (B0 + kk * B1 + kk * kk));
    const float gdsp = gvx * gblt;
    const float A1blt = ((2.0f * A0 - 2.0f * kk * kk) / (A0 + kk * A1 + kk * kk));
    const float A2blt = ((A0 - kk * A1 + kk * kk) / (A0 + kk * A1 + kk * kk));
    const float B1blt = ((2.0f * B0 - 2.0f * kk * kk) / (B0 + kk * B1 + kk * kk));
    const float B2blt = ((B0 - kk * B1 + kk * kk) / (B0 + kk * B1 + kk * kk));

    ks->st1.fs = ks->fs;
    ks->st1.a1 = -B1blt;  // negate because of filter implementation
    ks->st1.a2 = -B2blt;  // negate because of filter implementation
    ks->st1.b0 = 1.0f;
    ks->st1.b1 = A1blt;
    ks->st1.b2 = A2blt;
    ks->st1.gain = gdsp;

    // Stage 2 filter (order 1)
    // vovxblt = gdsp2*(X0blt*z1 + 1)/(Y0blt*z1 + 1)
    const float g2blt = (X0 + kk) / (Y0 + kk);
    const float gdsp2 = gvo * g2blt;
    const float X0blt = (X0 - kk) / (X0 + kk);
    const float Y0blt = (Y0 - kk) / (Y0 + kk);

    ks->st2.fs = ks->fs;
    ks->st2.a1 = -Y0blt;  // negate because of filter implementation
    ks->st2.a2 = 0.0f;
    ks->st2.b0 = 1.0f;
    ks->st2.b1 = X0blt;
    ks->st2.b2 = 0.0f;
    ks->st2.gain = gdsp2;
    //printf("b0= %f b1= %f b2= %f a1=%f a2= %f gain= %f\n", ks->st2.b0,ks->st2.b1, ks->st2.b2,ks->st2.a1,ks->st2.a2 ,ks->st2.gain );
    //printf("b0= %f b1= %f b2= %f a1=%f a2= %f gain= %f\n", ks->st1.b0,ks->st1.b1, ks->st1.b2,ks->st1.a1,ks->st1.a2 ,ks->st1.gain );

}

void kotstack_init(kot_stack *ks, float fs_) {
    // Sampling parameters
    ks->fs = fs_;

    // Powers of 10
    const float k = 1000.0f;
    const float n = 1.0e-9f;

    // Tonestack components
    ks->ri = 1.0f * k;
    ks->rc = 6.8f * k;
    ks->ro = 100.0f * k * 1000.0f * k / (100.0f * k + 1000.0f * k);
    ks->c1 = 10.0f * n;
    ks->c2 = 10.0f * n;
    ks->rtone = 25.0f * k;     // Tone Pot
    ks->rboost = 50.0f * k;      // Boost pot

    // Tone settings
    ks->tone_pot_pos = 0.2f;
    ks->boost_pot_pos = 0.25f;

    kotstack_compute_filter_coeffs(ks);

    ks->st1.x1 = 0.0f;
    ks->st1.x2 = 0.0f;
    ks->st1.y1 = 0.0f;
    ks->st1.y2 = 0.0f;

    ks->st2.x1 = 0.0f;
    ks->st2.x2 = 0.0f;
    ks->st2.y1 = 0.0f;
    ks->st2.y2 = 0.0f;

}

// User functions
void kotstack_set_tone(kot_stack *ks, float tone) {
    if (tone > 1.0f) {
        ks->tone_pot_pos = 0.0f;
    } else if (tone < 0.0f) {
        ks->tone_pot_pos = 1.0f;
    } else {
        ks->tone_pot_pos = 1.0f - tone;
    }

    kotstack_compute_filter_coeffs(ks);
}

void kotstack_set_boost(kot_stack *ks, float boost) {
    if (boost > 1.0f) {
        ks->boost_pot_pos = 1.0f;
    } else if (boost < 0.0f) {
        ks->boost_pot_pos = 0.0f;
    } else {
        ks->boost_pot_pos = boost;
    }

    kotstack_compute_filter_coeffs(ks);
}
