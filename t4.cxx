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

  SMatS a, b;
  SMat  c;

  std::default_random_engine      gen(0xbeef0133);
  std::normal_distribution<float> dis(1.0, 0.05);

  long64 count = 0;

init:
  for (int i = 0; i < M; ++i)
  {
    for (int j = i; j < M; ++j)
    {
      a(i,j) = dis(gen);
      b(i,j) = dis(gen);
    }
  }

  mpt.MPlexSym(0, 0).Assign(0, a.Array());
  mpt.MPlexSym(1, 0).Assign(0, b.Array());

  c = a * b;

  mpt.mult2_sym(1);

  if (++count % 100000 == 0)
    printf("Count = %lld\n", count);

  bool dump = false;

  for (int j = 0; j < M; ++j)
    for (int k = 0; k < M; ++k)
    {
      // There are occasional diffs up to 4.768372e-07 on host, very very
      // rarely on MIC. Apparently this is a rounding difference between AVX
      // and normal maths. On MIC it might be usage of FMA?
      if (std::abs(c(j,k) - mpt.MPlex(0, 0).At(0,j, k)) > 3e-7)
      {
        dump = true;
        printf("%d,%d d=%e (count = %lld)\n", j, k, c(j,k) - mpt.MPlex(0, 0).At(0, j, k), count);
      }
    }

  if (dump)
  {
    for (int i = 0; i < M; ++i)
    {
      for (int j = 0; j < M; ++j)
        printf("%8f ", c(i,j));
      printf("\n");
    }
    printf("\n");

    for (int i = 0; i < M; ++i)
    {
      for (int j = 0; j < M; ++j)
        printf("%8f ", mpt.MPlex(0, 0).At(0,i,j));
      printf("\n");
    }
    printf("\n");
  }

  goto init;

  return 0;
}
