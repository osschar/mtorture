#!/usr/bin/perl

use lib ".";

use GenMul;

### If you're going to run GMtest.cxx and you do some changes here
### you *MUST* bring DIM, DOM and pattern assumptions in sync!

my $DIM = 3;
my $DOM = 6;

$a = new GenMul::MatrixSym('name'=>'a', 'M'=>$DIM);
$a->set_pattern(<<"FNORD");
x
x 1 
x x x
FNORD

$b = new GenMul::Matrix('name'=>'b', 'M'=>$DIM, 'N'=>$DOM);
$b->set_pattern(<<"FNORD");
x x x x 0 x
x 1 x 1 0 x
x x x x 0 x
FNORD

$c = new GenMul::Matrix('name'=>'c', 'M'=>$DIM, 'N'=>$DOM);

$bt = new GenMul::MatrixTranspose($b);

$bt->print_info();
$bt->print_pattern();

# ----------------------------------------------------------------------

$m = new GenMul::Multiply;


open STD, ">multify.ah";
select STD;

$m->multiply_standard($a, $b, $c);

close STD;

# print "\n", '-' x 80, "\n\n";

open INT, ">multify_intr.ah";
select INT;

$m->multiply_intrinsic($a, $b, $c);

close INT;
