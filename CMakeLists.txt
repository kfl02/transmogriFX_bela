cmake_minimum_required(VERSION 3.14)
project(transmogriFX_bela)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)

add_executable(transmogriFX_bela
        src/adsr.cpp
        src/adsr.h
        src/biquad.cpp
        src/biquad.h
        src/envelope_filter.cpp
        src/envelope_filter.h
        src/eq.cpp
        src/eq.h
        src/fb_compressor.cpp
        src/fb_compressor.h
        src/flange.cpp
        src/flange.h
        src/iir_1pole.cpp
        src/iir_1pole.h
        src/inductorwah.cpp
        src/inductorwah.h
        src/klingon.cpp
        src/klingon.h
        src/kot_tonestack.cpp
        src/kot_tonestack.h
        src/lfo.cpp
        src/lfo.h
        src/overdrive.cpp
        src/overdrive.h
        src/pareq.cpp
        src/pareq.h
        src/phaser.cpp
        src/phaser.h
        src/read_vi_trace.cpp
        src/read_vi_trace.h
        src/reverb.cpp
        src/reverb.h
        src/sample_hold_modulator.cpp
        src/sample_hold_modulator.h
        src/Sustainer.cpp
        src/Sustainer.h
        src/svf.cpp
        src/svf.h
        src/trem.cpp
        src/trem.h
        src/constants.h)
