#!/usr/bin/perl

use lib ".";

use GenMul;
use warnings;

# Shared multiply object
$m = new GenMul::Multiply;

#-----------------------------------------------------------
# 3x3
$A = new GenMul::Matrix('name'=>'a', 'M'=>3, 'N'=>3);
$B = new GenMul::Matrix('name'=>'b', 'M'=>3, 'N'=>3);
$C = new GenMul::Matrix('name'=>'c', 'M'=>3, 'N'=>3);

$m->dump_multiply_std_and_intrinsic("Matrix33x33.ah",
                                    $A, $B, $C);

#-----------------------------------------------------------
# 3x3 sym
$A = new GenMul::MatrixSym('name'=>'a', 'M'=>3, 'N'=>3);
$B = new GenMul::MatrixSym('name'=>'b', 'M'=>3, 'N'=>3);
$C = new GenMul::Matrix   ('name'=>'c', 'M'=>3, 'N'=>3);

$m->dump_multiply_std_and_intrinsic("MatrixSym33x33.ah",
                                    $A, $B, $C);

#-----------------------------------------------------------
# 6x6
$A = new GenMul::Matrix('name'=>'a', 'M'=>6, 'N'=>6);
$B = new GenMul::Matrix('name'=>'b', 'M'=>6, 'N'=>6);
$C = new GenMul::Matrix('name'=>'c', 'M'=>6, 'N'=>6);

$m->dump_multiply_std_and_intrinsic("Matrix66x66.ah",
                                    $A, $B, $C);

#-----------------------------------------------------------
# 6x6 sym
$A = new GenMul::MatrixSym('name'=>'a', 'M'=>6, 'N'=>6);
$B = new GenMul::MatrixSym('name'=>'b', 'M'=>6, 'N'=>6);
$C = new GenMul::Matrix   ('name'=>'c', 'M'=>6, 'N'=>6);

$m->dump_multiply_std_and_intrinsic("MatrixSym66x66.ah",
                                    $A, $B, $C);
