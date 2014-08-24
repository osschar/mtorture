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



##############################
### updateParameters       ###
##############################

#declared first on its own because propErr sees many uses
my $propErr_M = 6;
$propErr = new GenMul::MatrixSym('name'=>'a', 'M'=>$propErr_M); #will have to remember to re'name' it based on location in function
$propErr->set_pattern(<<"FNORD");
x
x x
x x x
x x x x
x x x x x
x x x x x x
FNORD

my $propErrT_M = 6;
$propErrT = new GenMul::MatrixTranspose($propErr); #will have to remember to re'name' it based on location in function



### kalmanGain =  = propErr * (projMatrixT * resErrInv)
my $resErrInv_proj63_M = 6;
my $resErrInv_proj63_N = 3;
$resErrInv_proj63 = new GenMul::Matrix('name'=>'b', 'M'=>$resErrInv_proj63_M, 'N'=>$resErrInv_proj63_N);
# this pattern is just right that because of the zeroes you should be able to just feed in the 3x3 without first converting to a 6x3????
$resErrInv_proj63->set_pattern(<<"FNORD");
x x x
x x x
x x x
0 0 0
0 0 0
0 0 0
FNORD

$kalmanGain = new GenMul::Matrix('name'=>'c', 'M'=>$propErr_M, 'N'=>$resErrInv_proj63_N);

$m->dump_multiply_std_and_intrinsic("upParam_MultKalmanGain.ah",
									$propErr, $resErrInv_proj63, $kalmanGain);


### updatedErrs = propErr - propErr^T * simil * propErr
# Going to skip the subtraction for now
my $simil_M = 6;
$simil = new GenMul::MatrixSym('name'=>'a', 'M'=>$simil_M);
$simil->set_pattern(<<"FNORD");
x
x x
x x x
0 0 0 0
0 0 0 0 0
0 0 0 0 0 0
FNORD

$propErr->{name} = 'b';

my $temp_simil_x_propErr_M = 6;
my $temp_simil_x_propErr_N = 6;
$temp_simil_x_propErr = new GenMul::Matrix('name'=>'c', 'M'=>$temp_simil_x_propErr_M, 'N'=>$temp_simil_x_propErr_N);

$m->dump_multiply_std_and_intrinsic("upParam_simil_x_propErr.ah",
									 $simil, $propErr, $temp_simil_x_propErr);

$temp_simil_x_propErr->{name} = 'b';									 
$temp_simil_x_propErr->set_pattern(<<"FNORD");
x x x x x x
x x x x x x
x x x x x x
0 0 0 0 0 0
0 0 0 0 0 0
0 0 0 0 0 0
FNORD

#? This one is symmetric but the output can't handle it... need to fix
#$temp_propErrT_x_simil_propErr = new GenMul::MatrixSym('name'=>'c', 'M'=>$propErrT_M, 'N'=>$temp_simil_x_propErr_N);
$temp_propErrT_x_simil_propErr = new GenMul::MatrixSym('name'=>'c', 'M'=>$propErrT_M);

$m->dump_multiply_std_and_intrinsic("upParam_propErrT_x_simil_propErr.ah",
									$propErrT, $temp_simil_x_propErr, $temp_propErrT_x_simil_propErr);
									

$kalmanGain66 = new GenMul::Matrix('name'=>'a', 'M'=>$propErr_M, 'N'=>$propErr_M);
$kalmanGain66->set_pattern(<<"FNORD");
x x x 0 0 0
x x x 0 0 0
x x x 0 0 0
x x x 0 0 0
x x x 0 0 0
x x x 0 0 0
FNORD

$temp_kalmanGain_x_propErr = new GenMul::MatrixSym('name'=>'c', 'M'=>$propErrT_M);

$m->dump_multiply_std_and_intrinsic("upParam_kalmanGain_x_propErr.ah",
									$kalmanGain66, $propErr, $temp_kalmanGain_x_propErr);
