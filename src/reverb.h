#pragma once

// -----------------------------------------------------------------------
//
//  Copyright (C) 2003-2011 Fons Adriaensen <fons@linuxaudio.org>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// -----------------------------------------------------------------------


#ifndef REVERB_H
#define REVERB_H

#include "pareq.h"


// -----------------------------------------------------------------------


class Diff1 {
private:

    friend class Reverb;

    Diff1();

    ~Diff1();

    void init(int size, float c);

    void fini();

    float process(float x) {
        float z = _line[_i];
        x -= _c * z;
        _line[_i] = x;
        if (++_i == _size) _i = 0;
        return z + _c * x;
    }

    int _i;
    float _c;
    int _size;
    float *_line;
};


// -----------------------------------------------------------------------


class Filt1 {
private:

    friend class Reverb;

    Filt1() : _slo(0), _shi(0) {}

    ~Filt1() {}

    void set_params(float del, float tmf, float tlo, float wlo, float thi, float chi);

    float process(float x) {
        _slo += _wlo * (x - _slo) + 1e-10f;
        x += _glo * _slo;
        _shi += _whi * (x - _shi);
        return _gmf * _shi;
    }

    float _gmf;
    float _glo;
    float _wlo;
    float _whi;
    float _slo;
    float _shi;
};


// -----------------------------------------------------------------------


class Delay {
private:

    friend class Reverb;

    Delay();

    ~Delay();

    void init(int size);

    void fini();

    float read() {
        return _line[_i];
    }

    void write(float x) {
        _line[_i++] = x;
        if (_i == _size) _i = 0;
    }

    int _i;
    int _size;
    float *_line;
};


// -----------------------------------------------------------------------


class Vdelay {
private:

    friend class Reverb;

    Vdelay();

    ~Vdelay();

    void init(int size);

    void fini();

    void set_delay(int del);

    float read() {
        float x = _line[_ir++];
        if (_ir == _size) _ir = 0;
        return x;
    }

    void write(float x) {
        _line[_iw++] = x;
        if (_iw == _size) _iw = 0;
    }

    int _ir;
    int _iw;
    int _size;
    float *_line;
};


// -----------------------------------------------------------------------


class Reverb {
public:

    Reverb();

    ~Reverb();

    void init(float fsamp, bool ambis, size_t frame_size);

    void fini();

    void prepare(int n);

    void process(int n, float *inp[], float *out[]);

    void tick_mono(int frames, float *audio);

    void set_delay(float v) {
        _ipdel = v;
        _cntA1++;
    }

    void set_xover(float v) {
        _xover = v;
        _cntB1++;
    }

    void set_rtlow(float v) {
        _rtlow = v;
        _cntB1++;
    }

    void set_rtmid(float v) {
        _rtmid = v;
        _cntB1++;
        _cntC1++;
    }

    void set_fdamp(float v) {
        _fdamp = v;
        _cntB1++;
    }

    void set_opmix(float v) {
        _opmix = v;
        _cntC1++;
    }

    void set_rgxyz(float v) {
        _rgxyz = v;
        _cntC1++;
    }

    void set_eq1(float f, float g) { _pareq1.setparam(f, g); }

    void set_eq2(float f, float g) { _pareq2.setparam(f, g); }

    float *out_mono;

private:


    float _fsamp;
    bool _ambis;

    Vdelay _vdelay0;
    Vdelay _vdelay1;
    Diff1 _diff1[8];
    Filt1 _filt1[8];
    Delay _delay[8];

    volatile int _cntA1;
    volatile int _cntB1;
    volatile int _cntC1;
    int _cntA2;
    int _cntB2;
    int _cntC2;

    float _ipdel;
    float _xover;
    float _rtlow;
    float _rtmid;
    float _fdamp;
    float _opmix;
    float _rgxyz;

    float _g0, _d0;
    float _g1, _d1;

    Pareq _pareq1;
    Pareq _pareq2;

    static float _tdiff1[8];
    static float _tdelay[8];

    //tick_mono variables
    int _fragm;
    int _nsamp;
    float *inp[2];
    float *out[4];


};

// -----------------------------------------------------------------------

#endif
