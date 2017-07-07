#include "ArrayTest.h"

#include "common.h"

//==============================================================================

ArrayTest::ArrayTest(int n_array, int size)
{
  fA = new_aligned<float*>(n_array);
  fN = n_array;

  for (int i = 0; i < fN; ++i)
  {
    fA[i] = new_aligned<float>(size);
  }
}

ArrayTest::~ArrayTest()
{
  for (int i = 0; i < fN; ++i)
  {
    free_aligned(fA[i]);
  }
  free_aligned(fA);
}

//==============================================================================

Func_t ArrayTest::name_to_func(const std::string& name)
{
  // perl -ne 'if (m/^\s*long64\s+(\w[\w\d]*)\(int\s+\w+\);\s*$/) { print "  if (name.compare(\"$1\") == 0) return [&](int n) { return $1(n); };\n"; }' ArrayTest.h

  if (name.compare("copy") == 0) return [&](int n) { return copy(n); };
  if (name.compare("sum2") == 0) return [&](int n) { return sum2(n); };
  if (name.compare("sum2_sqr") == 0) return [&](int n) { return sum2_sqr(n); };
  if (name.compare("sum2_cube") == 0) return [&](int n) { return sum2_cube(n); };
  if (name.compare("sum2_quad") == 0) return [&](int n) { return sum2_quad(n); };
  if (name.compare("sum2_quint") == 0) return [&](int n) { return sum2_quint(n); };
  if (name.compare("sum3") == 0) return [&](int n) { return sum3(n); };
  if (name.compare("sum3_sqr") == 0) return [&](int n) { return sum3_sqr(n); };
  if (name.compare("sum3_cube") == 0) return [&](int n) { return sum3_cube(n); };
  if (name.compare("mul2") == 0) return [&](int n) { return mul2(n); };
  if (name.compare("mul3") == 0) return [&](int n) { return mul3(n); };
  if (name.compare("div2") == 0) return [&](int n) { return div2(n); };
  if (name.compare("div3") == 0) return [&](int n) { return div3(n); };
  if (name.compare("sum2_cube_sa") == 0) return [&](int n) { return sum2_cube_sa(n); };
  if (name.compare("sum2_quint_sa") == 0) return [&](int n) { return sum2_quint_sa(n); };
  if (name.compare("sum3_cube_sa") == 0) return [&](int n) { return sum3_cube_sa(n); };
  if (name.compare("sin2") == 0) return [&](int n) { return sin2(n); };
  if (name.compare("cos2") == 0) return [&](int n) { return cos2(n); };
  if (name.compare("sincos2") == 0) return [&](int n) { return sincos2(n); };
  if (name.compare("sincos2_tyl4") == 0) return [&](int n) { return sincos2_tyl4(n); };
  if (name.compare("sincos2_tyl6") == 0) return [&](int n) { return sincos2_tyl6(n); };
  if (name.compare("atan2") == 0) return [&](int n) { return atan2(n); };

  throw std::runtime_error("ArrayTest::name_to_func no match");
}

//==============================================================================

long64 ArrayTest::copy(int n)
{
  float *Z = fA[0];
  float *A = fA[1];

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME(n%16, 0);

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

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  ASSUME(n%16, 0);

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

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  ASSUME(n%16, 0);

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

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  ASSUME(n%16, 0);

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

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  ASSUME(n%16, 0);

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

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  ASSUME(n%16, 0);

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

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

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

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

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

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

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
// Just multiplication and division
//------------------------------------------------------------------------------

long64 ArrayTest::mul2(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  // float *C = fA[3];

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  // ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    Z[i] = A[i] * B[i];
  }

  return n;
}
long64 ArrayTest::mul3(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  float *C = fA[3];

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    Z[i] = A[i] * B[i] * C[i];
  }

  return 2 * n;
}

long64 ArrayTest::div2(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  // float *C = fA[3];

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  // ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    Z[i] = A[i] / B[i];
  }

  return n;
}

long64 ArrayTest::div3(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  float *C = fA[3];

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    Z[i] = (A[i] / B[i]) / C[i];
  }

  return 2 * n;
}

//------------------------------------------------------------------------------
// Same Array clones (store into A, not Z)
//------------------------------------------------------------------------------

long64 ArrayTest::sum2_cube_sa(int n)
{
  // float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];

  // ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  ASSUME(n%16, 0);

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

  // ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  ASSUME(n%16, 0);

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

  // ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

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

//------------------------------------------------------------------------------
// Trigonometric
//------------------------------------------------------------------------------

namespace
{
  inline void sincos4(const float x, float& sin, float& cos)
  {
    const float x2 = x*x;
    cos  = 1.f - 0.5f*x2 + 0.041666667f*x2*x2;
    sin  = x - 0.16666667f*x*x2;
  }

  inline void sincos6(const float x, float& sin, float& cos)
  {
    const float x2 = x*x;
    const float x3 = x*x2;
    cos  = 1.f - 0.5f*x2 + 0.041666667f*x2*x2 - 0.0013888889f*x3*x3;
    sin  = x - 0.16666667f*x3 + 0.0083333333f*x3*x2;
  }
}

long64 ArrayTest::sin2(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  // float *C = fA[3];

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  // ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    Z[i] = std::sin(A[i]) + std::sin(B[i]);
  }

  return 2 * n;
}

long64 ArrayTest::cos2(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  // float *C = fA[3];

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  // ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    Z[i] = std::cos(A[i]) + std::cos(B[i]);
  }

  return 2 * n;
}

long64 ArrayTest::sincos2(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  // float *C = fA[3];

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  // ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    Z[i] = std::sin(A[i]) + std::sin(B[i]) +
           std::cos(A[i]) + std::cos(B[i]);
  }

  return 4 * n;
}

long64 ArrayTest::sincos2_tyl4(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  // float *C = fA[3];

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  // ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    float sa, ca, sb, cb;
    sincos4(A[i], sa, ca);
    sincos4(B[i], sb, cb);
    Z[i] = sa + sb + ca + cb;
  }

  return 4 * n;
}

long64 ArrayTest::sincos2_tyl6(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  // float *C = fA[3];

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  // ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    float sa, ca, sb, cb;
    sincos6(A[i], sa, ca);
    sincos6(B[i], sb, cb);
    Z[i] = sa + sb + ca + cb;
  }

  return 4 * n;
}

long64 ArrayTest::atan2(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  // float *C = fA[3];

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  // ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    Z[i] = std::atan2(B[i], A[i]);
  }

  return 1 * n;
}

//------------------------------------------------------------------------------
// Test template
//------------------------------------------------------------------------------

/*
long64 ArrayTest::xxx(int n)
{
  float *Z = fA[0];
  float *A = fA[1];
  float *B = fA[2];
  // float *C = fA[3];

  ASSUME_ALIGNED(Z, 64);
  ASSUME_ALIGNED(A, 64);
  ASSUME_ALIGNED(B, 64);
  // ASSUME_ALIGNED(C, 64);
  ASSUME(n%16, 0);

#pragma simd
  for (int i = 0; i < n; ++i)
  {
    // Z[i] = A[i] + B[i];
  }

  return 666 * n;
}
*/
