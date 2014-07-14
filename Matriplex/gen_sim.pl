#!/usr/bin/perl

@Offs[3] = [ 0, 1, 3, 1, 2, 4, 3, 4, 5 ];
@Offs[4] = [ 0, 1, 3, 6, 1, 2, 4, 7, 3, 4, 5, 8, 6, 7, 8, 9 ];
@Offs[5] = [ 0, 1, 3, 6, 10, 1, 2, 4, 7, 11, 3, 4, 5, 8, 12, 6, 7, 8, 9, 13, 10, 11, 12, 13, 14 ];
@Offs[6] = [ 0, 1, 3, 6, 10, 15, 1, 2, 4, 7, 11, 16, 3, 4, 5, 8, 12, 17, 6, 7, 8, 9, 13, 18, 10, 11, 12, 13, 14, 19, 15, 16, 17, 18, 19, 20 ];

$PREF = "  ";
$BR   = " ";
$JOIN = "$BR";
$POST = " ";

$D   = 6;
@Off = @{$Offs[$D]};

# $E = $D;
$E = 3; # for Kalman gain multiply

$SYMMETRIC = 1;
$RAW = 0;

$IN  = "psErr";
$OUT = "outErr";

if ($SYMMETRIC)
{ 
  if ($RAW)
  {
    for (my $i = 0; $i < $D; ++$i)
    {
      for (my $j = 0; $j <= $i; ++$j)
      {
        my $x = $Off[$i * $D + $j];
        
        my $x1 = $Off[($i - 3) * $D + $j];
        my $x2 = $Off[$i * $D + ($j - 3)];

        my $xx = $Off[($i - 3) * $D + ($j - 3)];

        print "${PREF}B.fArray[$x * N + n] = A.fArray[$x * N + n]";
      
        if ($i < 3 && $j < 3)
        {
        }
        elsif ($i >= 3 && $j >= 3)
        {
          print " + p * (A.fArray[$x1 * N + n] + A.fArray[$x2 * N + n]) + psq * A.fArray[$xx * N + n]";
        }
        elsif ($i >= 3)
        {
          print " + p * A.fArray[$x1 * N + n]";
        }
        elsif ($j >= 3)
        {
          print " + p * A.fArray[$x2 * N + n]";
        }
        else
        {
          print "Go to hell ...\n";
          exit 1;
        }

        print ";\n";
      }
    }
  }
  else
  {
    for (my $i = 0; $i < $D; ++$i)
    {
      for (my $j = 0; $j <= $i; ++$j)
      {
        my $x = $Off[$i * $D + $j];
        
        my $x1 = $Off[($i - 3) * $D + $j];
        my $x2 = $Off[$i * $D + ($j - 3)];

        my $xx = $Off[($i - 3) * $D + ($j - 3)];

        print "${PREF}${IN}[$x * N + n] = A.fArray[$x * N + n]";
      
        if ($i < 3 && $j < 3)
        {
        }
        elsif ($i >= 3 && $j >= 3)
        {
          print " + p * (A.fArray[$x1 * N + n] + A.fArray[$x2 * N + n]) + psq * A.fArray[$xx * N + n]";
        }
        elsif ($i >= 3)
        {
          print " + p * A.fArray[$x1 * N + n]";
        }
        elsif ($j >= 3)
        {
          print " + p * A.fArray[$x2 * N + n]";
        }
        else
        {
          print "Go to hell ...\n";
          exit 1;
        }

        print ";\n";
      }
    }
  }
}
else
{
  print "Not implemented ...\n";
  exit 1;

  for (my $i = 0; $i < $D; ++$i)
  {
    for (my $j = 0; $j < $D; ++$j)
    {
      my $x = $i * $D + $j;
      print "${PREF}C.fArray[$x * N + n] =${POST}";

    }
  }

}
