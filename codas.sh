#!/usr/bin/env bash

# Quick test
export TEST_DURATION=0.25

# Precise test -- takes a long while
# export TEST_DURATION=5.0

########################################################################

run_min()
{
    echo Entering loop for tests with minimal number of operations per loop
    for test in sum2 sum3 mul2 mul3 div2 div3 ; do
        echo Running ${test}
        ./t1 $test > arr_${test}_${pfix}.rt
    done
}

run_sig()
{
    echo Entering loop for tests with significant number of operations per loop
    for test in sum2_sqr sum2_cube sum2_quad sum2_quint sum3_sqr sum3_cube ; do
        echo Running ${test}
        ./t1 $test > arr_${test}_${pfix}.rt
    done
}

run_trig()
{
    echo Entering loop for tests with trigonometric functions
    for test in sin2 sincos2 sincos2_tyl4 sincos2_tyl6 atan2 ; do
        echo Running ${test}
        ./t1 $test > arr_${test}_${pfix}.rt
    done
}

#-----------------------------------------------------------------------

# Example comparing compilation options:
run_vset()
{
    test=sum2_quad

    for vset in sse4.2 avx; do
        # Have to rebuild for every change
        make clean
        make OPT:="-O3 -m${vset}" t1
        echo Running ${test} with -m${vset}
        ./t1 ${test} > arr_${test}_${vset}.rt
    done
}

########################################################################
# Main part of the script
########################################################################

make distclean
make t1

pfix=O3

# Minimal computation tests
run_min;

# Significant computation test
run_sig;

# Trigonometric functions
# run_trig;

# Vectorization options
# run_vset;
