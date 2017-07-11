#!/usr/bin/env bash

# Alternative to editing Makefile, pass options to make.
# E.g. for my osx laptop:
# MAKE_DEFS="OPT=\"-O3\" CXX:=c++-mp-5"

make distclean
make $MAKE_DEFS t1

# Precise test -- takes a long while
#export TEST_DURATION=5.0

pfix=O3

echo Entering loop for tests with minimal number of operations per loop
for test in sum2 sum3 mul2 mul3 div2 div3 ; do
    echo Running ${test}
    ./t1 $test > arr_${test}_${pfix}.rt
done

echo Entering loop for tests with significant number of operations per loop
for test in sum2_sqr sum2_cube sum2_quad sum2_quint sum3_sqr sum3_cube ; do
    echo Running ${test}
    ./t1 $test > arr_${test}_${pfix}.rt
done

echo Entering loop for tests with trigonometric functions
for test in sin2 sincos2 sincos2_tyl4 sincos2_tyl6 atan2 ; do
    echo Running ${test}
    ./t1 $test > arr_${test}_${pfix}.rt
done

exit 0

# Example comparing compilation options:

test=sum2_quint

for vset in sse4.2 avx avx2; do
    # Have to rebuild for every change
    make clean
    make OPT:="-O3 -m${vset}" t1
    echo Running ${test} with -m${vset}
    ./t1 ${test} > arr_${test}_${vset}.rt
done
