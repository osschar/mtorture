#include "ArrayTest.h"

#include "common.h"

//==============================================================================

ArrayTest::ArrayTest(int n_array, int size)
{
  fA = new_sth<float*>(n_array);
  fN = n_array;

  for (int i = 0; i < fN; ++i)
  {
    fA[i] = new_sth<float>(size);
  }
}

ArrayTest::~ArrayTest()
{
  for (int i = 0; i < fN; ++i)
  {
    free_sth(fA[i]);
  }
  free_sth(fA);
}

long64 ArrayTest::copy(int n)
{
  float *Z = fA[0];
  float *A = fA[1];

  __assume_aligned(Z, 64);
  __assume_aligned(A, 64);
  __assume(n%16==0);

  for (int i = 0; i < n; ++i)
  {
    Z[i] = A[i];
  }

  return n;
}

long64 ArrayTest::sum2(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];

  __assume_aligned(Z, 64);
  __assume_aligned(A, 64);
  __assume_aligned(B, 64);
  __assume(n%16==0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    Z[i] = A[i] + B[i];
  }

  return n;
}

long64 ArrayTest::sum2_sqr(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];

  __assume_aligned(Z, 64);
  __assume_aligned(A, 64);
  __assume_aligned(B, 64);
  __assume(n%16==0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    Z[i] = A[i]*A[i] + 2*A[i]*B[i] + B[i]*B[i];
  }

  return 6 * n;
}

long64 ArrayTest::sum2_cube(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];

  __assume_aligned(Z, 64);
  __assume_aligned(A, 64);
  __assume_aligned(B, 64);
  __assume(n%16==0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    const float asqr = A[i]*A[i];
    const float bsqr = B[i]*B[i];

    Z[i] = A[i]*asqr + 3*(asqr*B[i] + A[i]*bsqr) + B[i]*bsqr;
  }

  return 10 * n;
}

long64 ArrayTest::sum2_quad(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];

  __assume_aligned(Z, 64);
  __assume_aligned(A, 64);
  __assume_aligned(B, 64);
  __assume(n%16==0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    const float asqr = A[i]*A[i], acub = A[i]*asqr;
    const float bsqr = B[i]*B[i], bcub = B[i]*bsqr;

    Z[i] = A[i]*acub + 4*(acub*B[i] + A[i]*bcub) + 6*asqr*bsqr + B[i]*bcub;
  }

  return 15 * n;
}

long64 ArrayTest::sum2_quint(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];

  __assume_aligned(Z, 64);
  __assume_aligned(A, 64);
  __assume_aligned(B, 64);
  __assume(n%16==0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    const float asqr = A[i]*A[i], acub = A[i]*asqr, aqud = asqr*asqr;
    const float bsqr = B[i]*B[i], bcub = B[i]*bsqr, bqud = bsqr*bsqr;

    Z[i] = A[i]*aqud + 5*(aqud*B[i] + A[i]*bqud) + 10*(acub*bsqr + asqr*bcub) + B[i]*bqud;
  }

  return 19 * n;
}

//------------------------------------------------------------------------------

long64 ArrayTest::sum3(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  float *C = fA[3];

  __assume_aligned(Z, 64);
  __assume_aligned(A, 64);
  __assume_aligned(B, 64);
  __assume_aligned(C, 64);
  __assume(n%16==0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    Z[i] = A[i] + B[i] + C[i];
  }

  return 2 * n;
}

long64 ArrayTest::sum3_sqr(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  float *C = fA[3];

  __assume_aligned(Z, 64);
  __assume_aligned(A, 64);
  __assume_aligned(B, 64);
  __assume_aligned(C, 64);
  __assume(n%16==0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    Z[i] = A[i]*A[i] + B[i]*B[i] + C[i]*C[i] +
        2*(A[i]*B[i] + A[i]*C[i] + C[i]*B[i]);
  }

  return 12 * n;
}

long64 ArrayTest::sum3_cube(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  float *C = fA[3];

  __assume_aligned(Z, 64);
  __assume_aligned(A, 64);
  __assume_aligned(B, 64);
  __assume_aligned(C, 64);
  __assume(n%16==0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    const float asqr = A[i]*A[i];
    const float bsqr = B[i]*B[i];
    const float csqr = C[i]*C[i];

    Z[i] = A[i]*asqr + B[i]*bsqr + C[i]*csqr + 6*A[i]*B[i]*C[i] +
      3*(asqr*(B[i] + C[i]) + bsqr*(A[i] + C[i]) + csqr*(A[i] + B[i]));
  }

  return 22 * n;
}

//------------------------------------------------------------------------------
// Same Array clones (store into A, not Z)
//------------------------------------------------------------------------------

long64 ArrayTest::sum2_cube_sa(int n)
{
  // float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];

  // __assume_aligned(Z, 64);
  __assume_aligned(A, 64);
  __assume_aligned(B, 64);
  __assume(n%16==0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    const float asqr = A[i]*A[i];
    const float bsqr = B[i]*B[i];

    A[i] = A[i]*asqr + 3*(asqr*B[i] + A[i]*bsqr) + B[i]*bsqr;
  }

  return 10 * n;
}

long64 ArrayTest::sum2_quint_sa(int n)
{
  // float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];

  // __assume_aligned(Z, 64);
  __assume_aligned(A, 64);
  __assume_aligned(B, 64);
  __assume(n%16==0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    const float asqr = A[i]*A[i], acub = A[i]*asqr, aqud = asqr*asqr;
    const float bsqr = B[i]*B[i], bcub = B[i]*bsqr, bqud = bsqr*bsqr;

    A[i] = A[i]*aqud + 5*(aqud*B[i] + A[i]*bqud) + 10*(acub*bsqr + asqr*bcub) + B[i]*bqud;
  }

  return 19 * n;
}

long64 ArrayTest::sum3_cube_sa(int n)
{
  // float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  float *C = fA[3];

  // __assume_aligned(Z, 64);
  __assume_aligned(A, 64);
  __assume_aligned(B, 64);
  __assume_aligned(C, 64);
  __assume(n%16==0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    const float asqr = A[i]*A[i];
    const float bsqr = B[i]*B[i];
    const float csqr = C[i]*C[i];

    A[i] = A[i]*asqr + B[i]*bsqr + C[i]*csqr + 6*A[i]*B[i]*C[i] +
      3*(asqr*(B[i] + C[i]) + bsqr*(A[i] + C[i]) + csqr*(A[i] + B[i]));
  }

  return 22 * n;
}
