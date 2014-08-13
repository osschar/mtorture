#!/usr/bin/perl

use lib ".";

use GenMul;

my $DIM = 3;
my $DOM = 6;

$a = new GenMul::MatrixSym('name'=>'a', 'M'=>$DIM);

$b = new GenMul::Matrix('name'=>'b', 'M'=>$DIM, 'N'=>$DOM);
$b->set_pattern(<<"FNORD");
1 0 0 x 0 0
0 1 0 0 x 0
0 0 1 0 0 x
FNORD

$c = new GenMul::Matrix('name'=>'c', 'M'=>$DIM, 'N'=>$DOM);

# ----------------------------------------------------------------------

$m = new GenMul::Multiply;

$m->multiply_standard($a, $b, $c);

# $m->multiply_intrinsic($a, $b, $c);
