#include "common.h"
#include "Timing.h"
#include "ArrayTest.h"

// See environment variable influences in common.h

#ifndef TEST_FUNC
#define TEST_FUNC sum2_cube
#endif

#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif

#ifndef THREAD_BINDING
#define THREAD_BINDING spread
#endif


int main()
{
  const int NT = NUM_THREADS;

  ArrayTest *ats[NT];

  for (int i = 0; i < NT; ++i)
  {
    ats[i] = new ArrayTest(4, g_n_vec_max);
  }

  Timing t([&](int n_vec)
           {
             return ats[0]->TEST_FUNC(n_vec);
           });


  t.print_tuple_header();
  // t.print_header();

  // Warm up
#pragma omp parallel for num_threads(NT), proc_bind(THREAD_BINDING)
  for (int i = 0; i < NT; ++i)
  {
    for (long64 l = 0; l < 4; ++l)
    {
      ats[i]->TEST_FUNC(g_n_vec_max);
    }
  }

  for (int n_vec = g_n_vec_min; n_vec <= g_n_vec_max; n_vec *= 2)
  {
    long64 n_loop = t.calibrate_loop(n_vec, g_test_duration);

    long64 n_ops = 0;

    t.start();

// WTH ... thread_binding ignored ?
// #pragma omp parallel for num_threads(NT), proc_bind(THREAD_BINDING), reduction(+:n_ops)
#pragma omp parallel for num_threads(NT), reduction(+:n_ops)
    for (int i = 0; i < NT; ++i)
    {
      for (long64 l = 0; l < n_loop; ++l)
      {
        n_ops += ats[i]->TEST_FUNC(n_vec);
      }
    }

    t.stop(n_ops);

    t.print(n_vec);
  }

  return 0;
}
