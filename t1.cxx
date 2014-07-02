#include "Timing.h"
#include "ArrayTest.h"

int main()
{
  const int n_vec_min  = 8;
  const int n_vec_max  = 64 * 1024 * 1024;

  ArrayTest at(4, n_vec_max);

  Timing t([&at](int n_vec)
           {
             return at.sum2(n_vec);
             // return at.sum2_cube(n_vec);
           });

  t.print_tuple_header();
  // t.print_header();

  // Warm up
  t.time_loop(n_vec_max, 4);

  for (int n_vec = n_vec_min; n_vec <= n_vec_max; n_vec *= 2)
  {

    t.auto_time_loop(n_vec, 1.0);

    t.print(n_vec);
  }

  return 0;
}
