#!/usr/bin/perl

use lib ".";
use Test;

########################################################################

my $t = new Test('Base'          => 't3',
                 'FmtHost'       => 'XXX_${test}_host.rt',
                 'FmtMic'        => 'XXX_${test}_mic.rt',

                 'RunHost'       => 1,
                 'RunMic'        => 1,

                 'MAKE_VARS'     => 'DEFS="-DMPT_DIM=3 -DMPT_SIZE=16"',

                 # Defaults for test duration and problem-size limits:
                 # 'TEST_DURATION' => 1,
                 # 'N_VEC_MIN'     => 1,
                 # 'N_VEC_MAX'     => 64 * 1024 * 1024,

                 # "Precise" test:
                 'TEST_DURATION' => 5,

                 # "Small" test setting for development
                 #'TEST_DURATION' => 0.5,
                 #'N_VEC_MIN'     => 1,
                 #'N_VEC_MAX'     => 16384,

                'EnvHost' => "",
                'EnvMic'  => "",
    );

@TESTS = qw(mult2_sym);
# Other tests:
# - mult2_3out mult2_3in            --- made no real difference
# - inv_cramer inv_cholesky         --- works only for 3x3
# - inv_cramer_sym inv_cholesky_sym --- works only for 3x3

for $dim (3)#, 6)
{
  for $size (8, 16, 32, 64)
  {
    $t->{'FmtHost'} = "\${test}_${dim}_${size}_host.rt";
    $t->{'FmtMic'}  = "\${test}_${dim}_${size}_mic.rt";

    $t->{MAKE_VARS} = "DEFS=\"-DMPT_DIM=$dim -DMPT_SIZE=$size\"";

    $t->make_clean();

    for $test (@TESTS)
    {
      print "Running test $test ...\n";

      $t->run_test($test);
    }
  }
}
