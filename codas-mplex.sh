#!/usr/bin/env bash

# Quick test
# export TEST_DURATION=0.25

# Precise test -- takes a long while
export TEST_DURATION=5.0

########################################################################

run_mplex_mults()
{
    pfx=mplx_mult2
    base_opts="-O3"

    #for dim in 6; do
    #    for size in 8; do

    for dim in 3 6
    do
        for size in 8 16 32
        do
            sfx=${dim}_${size}.rt
            opts="${base_opts} -DMPT_DIM=${dim} -DMPT_SIZE=${size}"
            
            echo XXXXX Runnung for OPTS ${opts}
            echo

            make clean
            make OPTS:="${opts}" t3

            ./t3 mult2 > ${pfx}_std_${sfx}
            ./t3 mult2_general > ${pfx}_general_${sfx}

            ./t3 mult2_sym > ${pfx}_sym_std_${sfx}
            ./t3 mult2_sym_general > ${pfx}_sym_general_${sfx}

            opts="${opts} -DMPLEX_USE_INTRINSICS"

            echo XXXXX Runnung with OPTS ${opts}
            echo

            make clean
            make OPTS:="${opts}" t3

            ./t3 mult2 > ${pfx}_intr_${sfx}
            ./t3 mult2_sym > ${pfx}_sym_intr_${sfx}
        done
    done
}

########################################################################
# Main part of the script
########################################################################

run_mplex_mults
