#!/usr/bin/perl

use lib "../Matriplex";

use GenMul;

my $DIM = 6;

### Propagate Helix To R -- final similarity, two ops.

# outErr = errProp * outErr * errPropT
#   outErr is symmetric


$errProp = new GenMul::Matrix('name'=>'a', 'M'=>$DIM, 'N'=>$DIM);
$errProp->set_pattern(<<"FNORD");
x x 0 x x 0
x x 0 x x 0
x x 1 x x x
x x 0 x x 0
x x 0 x x 0
0 0 0 0 0 1
FNORD

$outErr = new GenMul::MatrixSym('name'=>'b', 'M'=>$DIM, 'N'=>$DIM);

$temp   = new GenMul::Matrix('name'=>'c', 'M'=>$DIM, 'N'=>$DIM);


$errPropT = new GenMul::MatrixTranspose($errProp);
$errPropT->print_info();
$errPropT->print_pattern();

# ----------------------------------------------------------------------

$m = new GenMul::Multiply;

# outErr and c are just templates ...

$m->dump_multiply_std_and_intrinsic("MultHelixProp.ah",
                                    $errProp, $outErr, $temp);

$temp  ->{name} = 'b';
$outErr->{name} = 'c';

### XXX fix this ... in accordance with what is in Propagation.cc
$m->dump_multiply_std_and_intrinsic("MultHelixPropTransp.ah",
                                    $temp, $errPropT, $outErr);
