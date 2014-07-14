#!/usr/bin/perl

use lib ".";
use Test;

########################################################################

my $t = new Test('Base'          => 't1',
                 'FmtHost'       => 'arr_${test}_host.rt',
                 'FmtMic'        => 'arr_${test}_mic.rt',

                 'RunHost'       => 1,
                 'RunMic'        => 1,

                 # 'MAKE_VARS'     => 'OPT=-O3',

                 # Defaults for test duration and problem-size limits:
                 # 'TEST_DURATION' => 1,
                 # 'N_VEC_MIN'     => 1,
                 # 'N_VEC_MAX'     => 64 * 1024 * 1024,

                 # "Precise" test:
                 'TEST_DURATION' => 5,

                 # "Small" test setting for development
                 # 'TEST_DURATION' => 0.2,
                 # 'N_VEC_MIN'     => 1024,
                 # 'N_VEC_MAX'     => 8 * 1024,
    );

@TESTS = qw(sum2 sum2_sqr sum2_cube sum2_quad sum2_quint
            sum3 sum3_sqr sum3_cube);

$t->make_clean();

for $test (@TESTS)
{
  print "Running test $test ...\n";

  $t->run_test($test);
}
