#!/bin/bash
cd "$(dirname "$0")"
# compile using the recommended settings
gcc -Wall -std=c99 -msse2 -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -m64 -lm \
    tvb.c -o tvb
# create output dir (otherwise the program will crash)
mkdir "./output"
# # run simulation for test subject
# ./tvb param_set_1 ./test_input
# g++ -o check_fit check_fit.cpp ../gsl_build/lib/libgsl.a ../gsl_build/lib/libgslcblas.a -I/data/project/ei_development/tools/gsl_build/include -lm && ./check_fit