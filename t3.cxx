#include "common.h"
#include "Timing.h"
#include "MPlexTest.h"

// See environment variable influences in common.h

#ifndef TEST_FUNC
#define TEST_FUNC mult2
#endif

const int mp_vec_max = g_n_vec_max / (MPT_DIM * MPT_DIM * MPT_SIZE);

int main(int argc, char *argv[])
{
  MPlexTest mpt(3, 2, mp_vec_max);

  Func_t foo;

  if (argc == 1)
  {
    foo = [&mpt](int n_vec) { return mpt.TEST_FUNC(n_vec); };
  }
  else
  {
    foo = mpt.name_to_func(argv[1]);
  }
  
  Timing t(foo);

  t.print_tuple_header();
  // t.print_header();

  // Warm up
  t.time_loop(mp_vec_max, 4);

  for (int n_vec = g_n_vec_min; n_vec <= mp_vec_max; n_vec *= 2)
  {
    t.auto_time_loop(n_vec, g_test_duration);

    t.print(n_vec);
  }

  return 0;
}
