#!/usr/bin/env bash

make t1

# Limit benchmark extent / duration

#export N_VEC_MIN=32
#export N_VEC_MAX=65536
#export TEST_DURATION=0.5

# Precise test -- takes a long while
#export TEST_DURATION=5.0

echo Entering loop for tests with minimal number of operations per loop
for test in sum2 sum3 mul2 mul3 div2 div3 ; do
    echo Running ${test}
    ./t1 $test > arr_${test}_min.rt
done

echo Entering loop for tests with significant number of operations per loop
for test in sum2_sqr sum2_cube sum2_quad sum2_quint sum3_sqr sum3_cube ; do
    echo Running ${test}
    ./t1 $test > arr_${test}_sig.rt
done

echo Entering loop for tests with trigonometric functions
for test in sin2 sincos2 sincos2_tyl4 sincos2_tyl6 atan2 ; do
    echo Running ${test}
    ./t1 $test > arr_${test}_trig.rt
done
