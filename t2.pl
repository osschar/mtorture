#!/usr/bin/perl

use lib ".";
use Test;

########################################################################

my $t = new Test('Base'          => 't2',
                 'FmtHost'       => 'mt_${test}_host.rt',
                 'FmtMic'        => 'mt_${test}_mic.rt',

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
                 # 'TEST_DURATION' => 0.5,
                 # 'N_VEC_MIN'     => 1024,
                 # 'N_VEC_MAX'     => 128 * 1024,

                'EnvHost' => "KMP_AFFINITY=scatter",
                'EnvMic'  => "KMP_AFFINITY='verbose,granularity=fine,compact' KMP_PLACE_THREADS=1T",
    );

@TESTS = qw(sum2_cube sum2_quint sum3_cube);

$t->make_clean();

for $test (@TESTS)
{
  print "Running test $test ...\n";

  $t->run_test($test);
}
