#include "Timing.h"

#include <cstdio>

#ifdef __MIC__
long64 Timing::s_cpu_freq       = 1238094000 / 2;
int    Timing::s_vec_unit_width = 16;
#else
long64 Timing::s_cpu_freq       = 2000003000;
int    Timing::s_vec_unit_width = 8;
#endif

double Timing::ops_per_tick()
{
  return ((double) m_n_ops) / m_diff / s_cpu_freq;
}

double Timing::vector_utilization()
{
  return ops_per_tick() / s_vec_unit_width;
}

//------------------------------------------------------------------------------

void Timing::print_tuple_header()
{
  printf("%s/I:%s/D:%s:%s:%s:%s\n",
         "NVec", "Time", "Gops", "Gflops", "OpT", "VecUt");
}

void Timing::print_header()
{
  printf("%8s   %10s   %10s   %10s   %8s   %6s\n",
         "NVec", "Time", "Gops", "Gflops", "OpT", "VecUt");
}

void Timing::print(int n_vec)
{
  printf("%8d   %10.6f   %10.6f   %10.6f   %8.4f   %6.4f\n",
         n_vec, m_diff, Gops(), Gflops(), ops_per_tick(), vector_utilization());

}

//------------------------------------------------------------------------------

void Timing::time_loop(int n_vec, int n_loop)
{
  m_n_ops = 0;

  start();

  for (int i = 0; i< n_loop; ++i)
  {
    m_n_ops += m_func(n_vec);
  }

  stop();
}

void Timing::auto_time_loop(int n_vec, double run_time)
{
  long64 n_ops = m_func(n_vec);

  // Estimate run time at full efficiency, run for 5% of that.

  long64 n_loops = s_cpu_freq * s_vec_unit_width  / n_ops * (0.05 * run_time) + 1;

  // printf("Estimate n_loops = %lld\n", n_loops);

  time_loop(n_vec, n_loops);

  // Correct so that run-time will be ok.

  n_loops = run_time / m_diff * n_loops + 1;

  // printf("Run-time was %f, new estimate n_loops = %lld\n", m_diff, n_loops);

  time_loop(n_vec, n_loops);
}
