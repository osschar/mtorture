#!/usr/bin/perl

use lib ".";
use Test;

# for (my $nv = 1; $nv <= 64; $nv *= 2)
for $nh (50, 100)
{
for $nv (1, 16)
{
  Test::system_or_die("make clean");

  Test::system_or_die("make DEFS=\"-DMAX_HITS=$nh -DMPT_SIZE=$nv\" mkFit -j 24");

  system "echo -n \"NVEC=$nv NHIT=$nh \" >> vec-scan-host.txt ";
  Test::system_or_die("mkFit/mkFit >> vec-scan-host.txt");

  system "echo -n \"NVEC=$nv NHIT=$nh \" >> vec-scan-mic.txt ";
  Test::system_or_die("ssh mic0 ./mkFit-mic >> vec-scan-mic.txt");
}
}

# for $nt (1, 2, 4, 8, 12, 16, 20, 24, 30, 45, 60, 75, 90, 105, 120)
# for $nt (10)
# {
#   Test::system_or_die("make clean");

#   Test::system_or_die("make DEFS=\"-DNUM_THREADS=$nt\" mkFit -j 24");

#   system "echo -n \"NTHR=$nt \" >> thr1-scan-host.txt ";
#   system "echo -n \"NTHR=$nt \" >> thr1-scan-mic.txt ";
#   system "echo -n \"NTHR=$nt \" >> thr2-scan-host.txt ";
#   system "echo -n \"NTHR=$nt \" >> thr2-scan-mic.txt ";

#   Test::system_or_die("KMP_AFFINITY=scatter KMP_PLACE_THREADS=1T mkFit/mkFit >> thr1-scan-host.txt");
#   Test::system_or_die("KMP_AFFINITY=compact KMP_PLACE_THREADS=2T mkFit/mkFit >> thr2-scan-host.txt");

#   Test::system_or_die("ssh mic0 KMP_AFFINITY=scatter KMP_PLACE_THREADS=1T ./mkFit-mic >> thr1-scan-mic.txt");
#   Test::system_or_die("ssh mic0 KMP_AFFINITY=compact KMP_PLACE_THREADS=2T ./mkFit-mic >> thr2-scan-mic.txt");
# }

# for $nt (1, 2, 4, 8, 12, 16, 20, 24, 30, 45, 60, 75, 90, 105, 120)
# for $nt (1, 2, 4, 8, 10, 12, 14, 16, 18, 20, 22, 24)
# for $nt (20, 2)
# {
#   Test::system_or_die("make clean");

#   # Test::system_or_die("make DEFS=\"-DNUM_THREADS=$nt\" mkFit -j 24");
#   Test::system_or_die("make DEFS=\"-DNUM_THREADS=$nt -DMPT_SIZE=8\" mkFit -j 24");

#   for ($ii = 0; $ii < 20; ++$ii)
#   {

#     #system "echo -n \"NTHR=$nt VEC=8 \" >> thr1-scan-host.txt ";
#     system "echo -n \"NTHR=$nt VEC=8 \" >> thr2-scan-host.txt ";
#     #system "echo -n \"NTHR=$nt \" >> thr1-scan-mic.txt ";
#     #system "echo -n \"NTHR=$nt \" >> thr2-scan-mic.txt ";

#     #Test::system_or_die("ssh mic0 KMP_AFFINITY=scatter KMP_PLACE_THREADS=1T ./mkFit-mic >> thr1-scan-mic.txt");
#     #Test::system_or_die("ssh mic0 KMP_AFFINITY=compact KMP_PLACE_THREADS=2T ./mkFit-mic >> thr2-scan-mic.txt");
#     #Test::system_or_die("KMP_AFFINITY=scatter KMP_PLACE_THREADS=1T mkFit/mkFit >> thr1-scan-host.txt");
#     Test::system_or_die("KMP_AFFINITY=compact KMP_PLACE_THREADS=2T mkFit/mkFit >> thr2-scan-host.txt");
#   }
# }
