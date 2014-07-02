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

int ArrayTest::copy(int n)
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

int ArrayTest::sum2(int n)
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

int ArrayTest::sum2_sqr(int n)
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

int ArrayTest::sum2_cube(int n)
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
    Z[i] = A[i]*asqr + 3*asqr*B[i] + 3*A[i]*bsqr + B[i]*bsqr;
  }

  return 11 * n;
}

//------------------------------------------------------------------------------

int ArrayTest::sum3(int n)
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

  return n;
}

int ArrayTest::sum3_sqr(int n)
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
int ArrayTest::sum3_cube(int n)
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
