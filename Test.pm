package Test;

use Carp;

########################################################################
# Helper functions
########################################################################

sub exec_or_die
{
  # Execute command and return its stdout.

  my $cmd = shift;

  my $ret = `$cmd`;
  croak "$cmd died" if $?;
  return $ret
}

sub system_or_die
{
  # Execute command via "system" call.

  my $cmd = shift;

  system "$cmd" and croak "$cmd died";
}

sub spawn
{
  # Fork a new process, redirect output to $log. Retunrs pid of the forked
  # process.

  my $log = shift;
  # Rest of arguments passed to exec

  my $pid = fork();

  if ($pid == 0)
  {
    open STDOUT, '>', $log  or die "Can't redirect STDOUT: $!";
    # open STDERR, ">&STDOUT" or die "Can't dup STDOUT: $!";

    select STDERR; $| = 1; # make unbuffered
    select STDOUT; $| = 1; # make unbuffered

    setpriority(0, 0, $NICE) if $NICE;

    exec @_;
  }

  return $pid;
}

sub wait_all
{
  # Wait for all child processes to finish.

  {
    # print "Waiting now ... ";
    my $pid = wait();
    my $est = $pid != -1 ? $? >> 8 : 'n/a';
    # print "pid = $pid, status = $est\n";

    redo if $pid != -1;
  }
}

########################################################################
# Class Test
########################################################################

sub new
{
  my $proto = shift;
  my $class = ref($proto) || $proto;

  my $S = {@_};
  bless($S, $class);

  ### For sub classes:
  ### @ISA = ('Test');
  ### my $S = $proto->SUPER::new(@_);

  $S->{'ExeEnv'}  = join(' ', map { defined $S->{$_} ? "$_=$S->{$_}" : () }
                         qw(TEST_DURATION PRE_TEST_FRAC N_VEC_MIN N_VEC_MAX
                            KMP_AFFINITY KMP_PLACE_THREADS));

  $S->{'Src'}     = "$S->{'Base'}.cxx";
  $S->{'ExeHost'} = "$S->{'Base'}";
  $S->{'ExeMic'}  = "$S->{'Base'}-mic";

  my @all_exes;
  push @all_exes, $S->{'ExeHost'} if $S->{'RunHost'};
  push @all_exes, $S->{'ExeMic'}  if $S->{'RunMic'};
  $S->{'AllExes'} = join(" ", @all_exes);

  $S->{'CmdHost'} = "         $S->{'ExeEnv'} $S->{'EnvHost'} ./$S->{'ExeHost'}";
  $S->{'CmdMic'}  = "ssh mic0 $S->{'ExeEnv'} $S->{'EnvMic'}  ./$S->{'ExeMic'}";

  return $S;
}

sub print_vars
{
  my $S = shift;

  for $k (sort keys %$S)
  {
    printf("%12s   =>   %s\n", $k, $S->{$k});
  }
}

########################################################################

sub make_clean
{
  my $S = shift;

  system_or_die("make $S->{'MAKE_VARS'} clean");
}

sub run_test
{
  my ($S, $test) = @_;

  # Force recompilation of test executable with appropriate define of TEST_FUNC.
  system_or_die("make $S->{'MAKE_VARS'} USER_CPPFLAGS=\"-DTEST_FUNC=$test\" -W $S->{'Src'} $S->{'AllExes'}");

  my $ohost = eval "\"$S->{'FmtHost'}\"";
  my $omic  = eval "\"$S->{'FmtMic'}\"";

  spawn($ohost, $S->{'CmdHost'}) if $S->{'RunHost'};
  spawn($omic,  $S->{'CmdMic'})  if $S->{'RunMic'};

  wait_all();
}

########################################################################

1;
