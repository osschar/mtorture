#include "common.h"
#include "MPlexTest.h"

#include "Math/SMatrix.h"

// See environment variable influences in common.h

#ifndef TEST_FUNC
#define TEST_FUNC mult2
#endif

const int mp_len     = S ;
const int mp_vec_max = g_n_vec_max / (MPT_DIM * MPT_DIM * MPT_SIZE);

const int M = MPT_DIM;

typedef ROOT::Math::SMatrix<float, MPT_DIM>                                                  SMat;
typedef ROOT::Math::SMatrix<float, MPT_DIM, MPT_DIM, ROOT::Math::MatRepSym<float, MPT_DIM> > SMatS;

int main()
{
  MPlexTest mpt(3, 2, 1);

  SMatS a[16], b[16];
  SMat  c[16];

  std::default_random_engine      gen(0xbeef0133);
  std::normal_distribution<float> dis(1.0, 0.05);

  long64 count = 0;

init:
  for (int m = 0; m < 16; ++m)
  {
    for (int i = 0; i < M; ++i)
    {
      for (int j = i; j < M; ++j)
      {
        a[m](i,j) = dis(gen);
        b[m](i,j) = dis(gen);
      }
    }

    mpt.MPlexSym(0, 0).Assign(m, a[m].Array());
    mpt.MPlexSym(1, 0).Assign(m, b[m].Array());

    c[m] = a[m] * b[m];
  }

  mpt.mult2_sym(1);

  if (++count % 100000 == 0)
    printf("Count = %lld\n", count);

  for (int m = 0; m < 16; ++m)
  {
    bool dump = false;

    for (int j = 0; j < M; ++j)
    {
      for (int k = 0; k < M; ++k)
      {
        // There are occasional diffs up to 4.768372e-07 on host, very very
        // rarely on MIC. Apparently this is a rounding difference between AVX
        // and normal maths. On MIC it might be usage of FMA?
        // The above was for 3x3.
        // For 6x6 practically all elements differ by 4.768372e-07, some
        // by 9.536743e-07.
        if (std::abs(c[m](j,k) - mpt.MPlex(0, 0).At(m, j, k)) > 5e-7)
        {
          dump = true;
          printf("M=%d  %d,%d d=%e (count = %lld)\n", m, j, k, c[m](j,k) - mpt.MPlex(0, 0).At(m, j, k), count);
        }
      }
    }

    if (dump && false)
    {
      printf("\n");
      for (int i = 0; i < M; ++i)
      {
        for (int j = 0; j < M; ++j)
          printf("%8f ", c[m](i,j));
        printf("\n");
      }
      printf("\n");

      for (int i = 0; i < M; ++i)
      {
        for (int j = 0; j < M; ++j)
          printf("%8f ", mpt.MPlex(0, 0).At(m, i, j));
        printf("\n");
      }
      printf("\n");
    }
    if (dump)
    {
      printf("\n");
    }
  }

  goto init;

  return 0;
}
