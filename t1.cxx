#include "common.h"
#include "Timing.h"
#include "ArrayTest.h"

// See environment variable influences in common.h

#ifndef TEST_FUNC
#define TEST_FUNC sum2_cube
#endif

int main(int argc, char *argv[])
{
  ArrayTest at(4, g_n_vec_max);

  Func_t foo;

  if (argc == 1)
  {
    foo = [&at](int n_vec) { return at.TEST_FUNC(n_vec); };
  }
  else
  {
    foo = at.name_to_func(argv[1]);
  }
  
  Timing t(foo);
  
  t.print_tuple_header();
  // t.print_header();

  // Warm up
  t.time_loop(g_n_vec_max, 4);

  for (int n_vec = g_n_vec_min; n_vec <= g_n_vec_max; n_vec *= 2)
  {
    t.auto_time_loop(n_vec, g_test_duration);

    t.print(n_vec);
  }

  return 0;
}
