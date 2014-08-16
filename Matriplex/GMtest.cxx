#include "Math/SMatrix.h"

#include "MatriplexSym.h"

#include <random>

/*
# Generate .ah files (make sure DIM, DOM and pattern match):
  ./GMtest.pl
# Compile host:
  icc -std=gnu++11 -openmp -mavx -O3 -I.. GMtest.cxx -o GMtest
# Compile MIC:
  icc -std=gnu++11 -openmp -mmic -O3 -I.. GMtest.cxx -o GMtest-mic && scp GMtest-mic mic0:
*/

typedef long long long64;

const int N   = 16;

const int DIM =  3;
const int DOM =  6;

#ifdef MIC_INTRINSICS
const float CMP_EPS = 2e-7;
#else
const float CMP_EPS = 4e-7;
#endif

typedef ROOT::Math::SMatrix<float, DIM, DOM>                                     SMatX;
typedef ROOT::Math::SMatrix<float, DIM, DIM, ROOT::Math::MatRepSym<float, DIM> > SMatS;

typedef Matriplex::Matriplex   <float, DIM, DOM, N>   MPlexX;
typedef Matriplex::MatriplexSym<float, DIM,      N>   MPlexS;

void Multify(const MPlexS& A, const MPlexX& B, MPlexX& C)
{
  typedef float T;

   const T *a = A.fArray; __assume_aligned(a, 64);
   const T *b = B.fArray; __assume_aligned(b, 64);
         T *c = C.fArray; __assume_aligned(c, 64);

#ifdef MIC_INTRINSICS

   for (int n = 0; n < N; n += 64 / sizeof(T))
   {
#include "multify_intr.ah"
   }

#else

#pragma simd
   for (int n = 0; n < N; ++n)
   {
#include "multify.ah"
   }

#endif
}

int main()
{
  SMatS  a[N];
  SMatX  b[N], c[N];

  MPlexS A;
  MPlexX B, C;

  std::default_random_engine      gen(0xbeef0133);
  std::normal_distribution<float> dis(1.0, 0.05);

  long64 count = 1;

init:

  for (int m = 0; m < N; ++m)
  {
    for (int i = 0; i < 3; ++i)
    {
      for (int j = i; j < 6; ++j)
      {
        if (j < DIM)  a[m](i,j) = dis(gen);

        b[m](i,j) = dis(gen);
      }
    }

    // Enforce pattern from GMtest.pl
    a[m](1, 1) = 1;
    b[m](0, 4) = 0;
    b[m](1, 1) = 1;
    b[m](1, 3) = 1;
    b[m](1, 4) = 0;
    b[m](2, 4) = 0;

    A.CopyIn(m, a[m].Array());
    B.CopyIn(m, b[m].Array());

    c[m] = a[m] * b[m];
  }

  Multify(A, B, C);

  for (int m = 0; m < N; ++m)
  {
    bool dump = false;

    for (int j = 0; j < DIM; ++j)
    {
      for (int k = 0; k < DOM; ++k)
      {
        // There are occasional diffs up to 4.768372e-07 on host, very very
        // rarely on MIC. Apparently this is a rounding difference between AVX
        // and normal maths. On MIC it might be usage of FMA?
        // The above was for 3x3.
        // For 6x6 practically all elements differ by 4.768372e-07, some
        // by 9.536743e-07.
        if (std::abs(c[m](j,k) - C.At(m, j, k)) > CMP_EPS)
        {
          dump = true;
          printf("M=%d  %d,%d d=%e (count = %lld)\n", m, j, k, c[m](j,k) - C.At(m, j, k), count);
        }
      }
    }

    if (dump && false)
    {
      printf("\n");
      for (int i = 0; i < DIM; ++i)
      {
        for (int j = 0; j < DOM; ++j)
          printf("%8f ", c[m](i,j));
        printf("\n");
      }
      printf("\n");

      for (int i = 0; i < DIM; ++i)
      {
        for (int j = 0; j < DOM; ++j)
          printf("%8f ", C.At(m, i, j));
        printf("\n");
      }
      printf("\n");
    }
    if (dump)
    {
      printf("\n");
    }
  }

  ++count;
  goto init;

  return 0;
}
