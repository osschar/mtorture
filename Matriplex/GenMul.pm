########################################################################
########################################################################
# MATRIX CLASSES
########################################################################
########################################################################

########################################################################
# MBase -- matrix base class
########################################################################

package GenMul::MBase;

use Carp;

# Required arguments:
# - M
# - N (for non-symmetric matrices)
# - name: name of array
#
# Created members:
# - class
#
#
# Input matrix pattern can be set via function set_pattern(). The argument is
# a white-space separated string of x, 0, 1, describing the matrix elements.
# For symmetric matrices lower-left triangle must be given.
# Support for -1 could be added (but isn't trivial (unless unary - changes 
# the preceeding addition into subtraction; also this is a tough call
# for intrinsics)).
#
# Pattern could also be set for output matrix but is currently not supported.

sub new
{
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my $S = {@_};
  bless($S, $class);

  # M, N checked in concrete classes

  croak "name must be set" unless defined $S->{name};

  $S->{class} = $class;

  return $S;
}

sub set_pattern
{
  my ($S, $pstr) = @_;

  @{$S->{pattern}} = split /\s+/, $pstr;

  croak "set_pattern number of entries does not match matrix size"
      unless scalar @{$S->{pattern}} == $S->mat_size();

  croak "set_pattern() input string contains invalid entry"
      if grep {$_ !~ /0|1|x/} @{$S->{pattern}};
}

sub mat_size
{
  die "max_size() should be overriden in concrete matrix class";
}

sub idx
{
  die "idx() should be overriden in concrete matrix class";
}

sub reg_name
{
  my ($S, $idx) = @_;

  return "$S->{name}_${idx}";
}

sub print_info
{
  my ($S) = @_;

  print "Class='$S->{class}', M=$S->{M}, N=$S->{N}\n";
}

########################################################################
# Matrix -- standard MxN matrix
########################################################################

package GenMul::Matrix; @ISA = ('GenMul::MBase');

use Carp;

sub new
{
  my $proto = shift;
  my $S = $proto->SUPER::new(@_);

  croak "M not set for $S->{class}" unless defined $S->{M};

  croak "N not set for $S->{class}" unless defined $S->{N};

  return $S;
}

sub mat_size
{
  my ($S) = @_;

  return $S->{M} * $S->{N};
}

sub idx
{
  my ($S, $i, $j) = @_;

  die "$S->{class}::idx() i out of range"
      if $i < 0 or $i >= $S->{M};

  die "$S->{class}::idx() j out of range"
      if $j < 0 or $j >= $S->{N};

  return $i * $S->{N} + $j;
}

########################################################################
# MatrixSym -- symmetric square matrix
########################################################################

package GenMul::MatrixSym; @ISA = ('GenMul::MBase');

use Carp;

# Offsets converting from full matrix indices to symmetric ones:
my @Offs;
@Offs[3] = [ 0, 1, 3, 1, 2, 4, 3, 4, 5 ];
@Offs[4] = [ 0, 1, 3, 6, 1, 2, 4, 7, 3, 4, 5, 8, 6, 7, 8, 9 ];
@Offs[5] = [ 0, 1, 3, 6, 10, 1, 2, 4, 7, 11, 3, 4, 5, 8, 12, 6, 7, 8, 9, 13, 10, 11, 12, 13, 14 ];
@Offs[6] = [ 0, 1, 3, 6, 10, 15, 1, 2, 4, 7, 11, 16, 3, 4, 5, 8, 12, 17, 6, 7, 8, 9, 13, 18, 10, 11, 12, 13, 14, 19, 15, 16, 17, 18, 19, 20 ];

sub new
{
  my $proto = shift;
  my $S = $proto->SUPER::new(@_);

  croak "M not set for $S->{class}" unless defined $S->{M};

  croak "N should not be set or should be equal to M for $S->{class}"
      if defined $S->{N} and $S->{N} != $S->{M};

  die "Offset array not defined for this dimension"
      unless defined @Offs[$S->{M}];

  die "Offset array of wrong dimension"
      unless scalar @{$Offs[$S->{M}]} == $S->{M} * $S->{M};

  $S->{N} = $S->{M} unless defined $S->{N};

  return $S;
}

sub mat_size
{
  my ($S) = @_;

  return ($S->{M} + 1) * $S->{M} / 2;
}

sub idx
{
  my ($S, $i, $j) = @_;

  die "$S->{class}::idx() i out of range"
      if $i < 0 or $i >= $S->{M};

  die "$S->{class}::idx() j out of range"
      if $j < 0 or $j >= $S->{N};

  return $Offs[$S->{M}][$i * $S->{N} + $j];
}


########################################################################
########################################################################
# CODE GENERATION CLASSES
########################################################################
########################################################################

package GenMul::Multiply;

use Carp;
use Scalar::Util 'blessed';

use warnings;


sub new
{
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my $S = {@_};
  bless($S, $class);

  $S->{prefix}  = "      " unless defined $S->{prefix};
  $S->{vectype} = "__m512" unless defined $S->{vectype};

  $S->{class} = $class;

  return $S;
}

sub check_multiply_arguments
{
  my ($S, $a, $b, $c) = @_;

  croak "Input a is not a GenMul::MBase"
      unless blessed $a and $a->isa("GenMul::MBase");

  croak "Input b is not a GenMul::MBase"
      unless blessed $b and $b->isa("GenMul::MBase");

  croak "Input c is not a GenMul::MBase"
      unless blessed $c and $c->isa("GenMul::MBase");

  croak "Input matrices a and b not compatible"
      unless $a->{N} == $b->{M};

  croak "Result matrix c of wrong dimensions"
      unless $c->{M} == $a->{M} and $c->{N} == $b->{N};

  croak "Result matrix c should not be symmetric (or implement this case in GenMul code)"
      if $c->isa("GenMul::MatrixSym");
}

sub push_out
{
  my $S = shift;

  push @{$S->{out}}, join "", @_;
}

sub unshift_out
{
  my $S = shift;

  unshift @{$S->{out}}, join "", @_;
}

sub handle_all_zeros_ones
{
  my ($S, $zeros, $ones) = @_;

  if ($zeros or $ones)
  {
    $S->unshift_out("");

    $S->unshift_out("$S->{vectype} all_ones  = { ", join(", ", (1) x 16), " };")
        if $ones;

    $S->unshift_out("$S->{vectype} all_zeros = { ", join(", ", (0) x 16), " };")
        if $zeros;
  }
}

# ----------------------------------------------------------------------

sub generate_addend_standard
{
  my ($S, $a, $aidx, $b, $bidx) = @_;

  my $apat = defined $a->{pattern} ? $a->{pattern}[$aidx] : 'x';
  my $bpat = defined $b->{pattern} ? $b->{pattern}[$bidx] : 'x';

  return undef if $apat eq '0' or  $bpat eq '0';
  return "1"   if $apat eq '1' and $bpat eq '1';

  my $astr = sprintf "$a->{name}\[%2d*N+n]", $aidx;
  my $bstr = sprintf "$b->{name}\[%2d*N+n]", $bidx;

  return $astr if $bpat eq '1';
  return $bstr if $apat eq '1';

  return "${astr}*${bstr}";
}

sub multiply_standard
{
  # Standard mutiplication - outputs unrolled C code, one line
  # per target matrix element.
  # Arguments: a, b, c   -- all GenMul::MBase with right dimensions.
  # Does:      c = a * b

  check_multiply_arguments(@_);

  my ($S, $a, $b, $c) = @_;

  for (my $i = 0; $i < $a->{M}; ++$i)
  {
    for (my $j = 0; $j < $b->{N}; ++$j)
    {
      my $x = $c->idx($i, $j);

      printf "$S->{prefix}$c->{name}\[%2d*N+n\] = ", $x;

      my @sum;

      for (my $k = 0; $k < $a->{N}; ++$k)
      {
        my $iko = $a->idx($i, $k);
        my $kjo = $b->idx($k, $j);

        my $addend = $S->generate_addend_standard($a, $iko, $b, $kjo);

        push @sum, $addend if defined $addend;
      }
      if (@sum)
      {
        print join(" + ", @sum), ";";
      }
      else
      {
        print "0;"
      }
      print "\n";
    }
  }
}

# ----------------------------------------------------------------------

sub load_if_needed
{
  my ($S, $mat, $idx, $arc) = @_;

  my $reg = $mat->reg_name(${idx});

  if ($arc->[$idx] == 0)
  {
    $S->push_out("$S->{vectype} ${reg} = LD($mat->{name}, $idx);");
    ++$S->{tick};
  }

  ++$arc->[$idx];

  return $reg;
}

sub store
{
  my ($S, $mat, $idx) = @_;

  my $reg = $mat->reg_name(${idx});

  $S->push_out("ST($mat->{name}, ${idx}, ${reg});");

  return $reg;
}

sub multiply_intrinsic
{
  check_multiply_arguments(@_);

  my ($S, $a, $b, $c) = @_;

  $S->{tick} = 0;

  $S->{out}  = [];

  # Counts of use. For a and b to fetch, for c to assign / add / mult / fma.
  # cc is used as tick at which store can be performed afterwards.
  my (@ac, @bc, @cc, @to_store);

  @ac = (0) x $a->mat_size();
  @bc = (0) x $b->mat_size();
  @cc = (0) x $c->mat_size();

  my $need_all_zeros = 0;
  my $need_all_ones  = 0;

  for (my $i = 0; $i < $a->{M}; ++$i)
  {
    for (my $k = 0; $k < $a->{N}; ++$k)
    {
      for (my $j = 0; $j < $b->{N}; ++$j)
      {
        my $x   = $c->idx($i, $j);
        my $iko = $a->idx($i, $k);
        my $kjo = $b->idx($k, $j);

        ### XXXX Need an intermediate check here for what operation to do.
        ### Can have:
        ### - add even for later operations
        ### - add 1 to all elements (but this should be rare).
        ### - assign on first operation (also with 1 as argument).
        ###
        ### ! All ones should be declared in the beginning, if needed.
        ### ! Dump all into a string, than prefix that if needed.

        my $apat = defined $a->{pattern} ? $a->{pattern}[$iko] : 'x';
        my $bpat = defined $b->{pattern} ? $b->{pattern}[$kjo] : 'x';

        if ($apat ne '0' and $bpat ne '0')
        {
          my ($areg, $breg, $sreg);

          if ($apat eq '1' and $bpat eq '1')
          {
            $need_all_ones = 1;
            $sreg = "all_ones";
          }
          elsif ($bpat eq '1')
          {
            $sreg = $S->load_if_needed($a, $iko, \@ac);
          }
          elsif ($apat eq '1')
          {
            $sreg = $S->load_if_needed($b, $kjo, \@bc);
          }
          else
          {
            $areg = $S->load_if_needed($a, $iko, \@ac);
            $breg = $S->load_if_needed($b, $kjo, \@bc);
          }

          my $creg = $c->reg_name($x);

          if ($cc[$x] == 0)
          {
            my $op = defined $sreg ? "${sreg}" : "MUL(${areg}, ${breg})";

            $S->push_out("$S->{vectype} ${creg} = ", $op, ";");
          }
          else
          {
            my $op = defined $sreg ?
                "ADD(${sreg}, ${creg})" :
                "FMA(${areg}, ${breg}, ${creg})";

            $S->push_out("${creg} = ", $op, ");");
          }

          ++$cc[$x];
          ++$S->{tick};
        }

        if ($k + 1 == $a->{N})
        {
          if ($cc[$x] == 0)
          {
            $need_all_zeros = 1;
            ### XXX emit store all_zeros into creg
          }
          else
          {
            $cc[$x] = $S->{tick} + 4; #### Will be ready to store in 4 cycles. Really 4?
            push @to_store, $x;
          }
        }

        # Try to store the finished ones.
        while (1)
        {
          last unless @to_store;
          my $s = $to_store[0];
          last if $S->{tick} < $cc[$s];

          $S->store($c, $s);
          shift @to_store;
          ++$S->{tick};
        }

      }

      $S->push_out("") unless $i + 1 == $a->{M} and $k + 1 == $a->{N};
    }
  }

  for my $s (@to_store)
  {
    $S->store($c, $s);

    ++$S->{tick};
  }

  $S->handle_all_zeros_ones($need_all_zeros, $need_all_ones);

  for (@{$S->{out}})
  {
    print $S->{prefix} unless /^$/;
    print;
    print "\n";
  }
}

########################################################################
########################################################################
# THE END
########################################################################
########################################################################

1;
