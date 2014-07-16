#ifndef Matriplex_H
#define Matriplex_H

#include "MatriplexCommon.h"

namespace Matriplex
{

//------------------------------------------------------------------------------

template<typename T, idx_t D1, idx_t D2, idx_t N>
class Matriplex
{
public:
   typedef T value_type;

   enum
   {
      /// return no. of matrix rows
      kRows = D1,
      /// return no. of matrix columns
      kCols = D2,
      /// return no of elements: rows*columns
      kSize = D1 * D2,
      /// size of the whole matriplex
      kTotSize = N * kSize
   };

   T fArray[kTotSize] __attribute__((aligned(64)));

   Matriplex() {}
   Matriplex(T v) { SetVal(v); }

   void SetVal(T v)
   {
      for (idx_t i = 0; i < kTotSize; ++i)
      {
         fArray[i] = v;
      }
   }

   T& At(idx_t n, idx_t i, idx_t j) { return fArray[(i * D2 + j) * N + n]; }

   T& operator()(idx_t n, idx_t i, idx_t j) { return fArray[(i * D2 + j) * N + n]; }

   void Assign(idx_t n, T *arr)
   {
      for (idx_t i = n; i < kTotSize; i += N)
      {
         fArray[i] = *(arr++);
      }
   }
};


template<typename T, idx_t D1, idx_t D2, idx_t N> using MPlex = Matriplex<T, D1, D2, N>;


//==============================================================================
// Multiplications
//==============================================================================

template<typename T, idx_t D1, idx_t D2, idx_t D3, idx_t N>
void Multiply(const MPlex<T, D1, D2, N>& A,
              const MPlex<T, D2, D3, N>& B,
              MPlex<T, D1, D3, N>& C)
{
   for (idx_t i = 0; i < D1; ++i)
   {
      for (idx_t j = 0; j < D3; ++j)
      {
         const idx_t ijo = N * (i * D3 + j);

         for (idx_t n = 0; n < N; ++n)
         {
            C.fArray[ijo + n] = 0;
         }

         //#pragma omp simd collapse(2)
         for (idx_t k = 0; k < D2; ++k)
         {
            const idx_t iko = N * (i * D2 + k);
            const idx_t kjo = N * (k * D3 + j);

#pragma simd
            for (idx_t n = 0; n < N; ++n)
            {
               // C.fArray[i, j, n] += A.fArray[i, k, n] * B.fArray[k, j, n];
               C.fArray[ijo + n] += A.fArray[iko + n] * B.fArray[kjo + n];
            }
         }
      }
   }
}

template<typename T, idx_t D1, idx_t D2, idx_t D3, idx_t N>
void MultiplyUnrolled(const MPlex<T, D1, D2, N>& A,
                      const MPlex<T, D2, D3, N>& B,
                      MPlex<T, D1, D3, N>& C)
{
   if (D1 == 3)
#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      C.fArray[ 0*N+n] = A.fArray[ 0*N+n]*B.fArray[ 0*N+n] + A.fArray[ 1*N+n]*B.fArray[ 3*N+n] + A.fArray[ 2*N+n]*B.fArray[ 6*N+n];
      C.fArray[ 1*N+n] = A.fArray[ 0*N+n]*B.fArray[ 1*N+n] + A.fArray[ 1*N+n]*B.fArray[ 4*N+n] + A.fArray[ 2*N+n]*B.fArray[ 7*N+n];
      C.fArray[ 2*N+n] = A.fArray[ 0*N+n]*B.fArray[ 2*N+n] + A.fArray[ 1*N+n]*B.fArray[ 5*N+n] + A.fArray[ 2*N+n]*B.fArray[ 8*N+n];
      C.fArray[ 3*N+n] = A.fArray[ 3*N+n]*B.fArray[ 0*N+n] + A.fArray[ 4*N+n]*B.fArray[ 3*N+n] + A.fArray[ 5*N+n]*B.fArray[ 6*N+n];
      C.fArray[ 4*N+n] = A.fArray[ 3*N+n]*B.fArray[ 1*N+n] + A.fArray[ 4*N+n]*B.fArray[ 4*N+n] + A.fArray[ 5*N+n]*B.fArray[ 7*N+n];
      C.fArray[ 5*N+n] = A.fArray[ 3*N+n]*B.fArray[ 2*N+n] + A.fArray[ 4*N+n]*B.fArray[ 5*N+n] + A.fArray[ 5*N+n]*B.fArray[ 8*N+n];
      C.fArray[ 6*N+n] = A.fArray[ 6*N+n]*B.fArray[ 0*N+n] + A.fArray[ 7*N+n]*B.fArray[ 3*N+n] + A.fArray[ 8*N+n]*B.fArray[ 6*N+n];
      C.fArray[ 7*N+n] = A.fArray[ 6*N+n]*B.fArray[ 1*N+n] + A.fArray[ 7*N+n]*B.fArray[ 4*N+n] + A.fArray[ 8*N+n]*B.fArray[ 7*N+n];
      C.fArray[ 8*N+n] = A.fArray[ 6*N+n]*B.fArray[ 2*N+n] + A.fArray[ 7*N+n]*B.fArray[ 5*N+n] + A.fArray[ 8*N+n]*B.fArray[ 8*N+n];
   }
   else if (D1 == 6)
#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      C.fArray[ 0*N+n] = A.fArray[ 0*N+n]*B.fArray[ 0*N+n] + A.fArray[ 1*N+n]*B.fArray[ 6*N+n] + A.fArray[ 2*N+n]*B.fArray[12*N+n] + A.fArray[ 3*N+n]*B.fArray[18*N+n] + A.fArray[ 4*N+n]*B.fArray[24*N+n] + A.fArray[ 5*N+n]*B.fArray[30*N+n];
      C.fArray[ 1*N+n] = A.fArray[ 0*N+n]*B.fArray[ 1*N+n] + A.fArray[ 1*N+n]*B.fArray[ 7*N+n] + A.fArray[ 2*N+n]*B.fArray[13*N+n] + A.fArray[ 3*N+n]*B.fArray[19*N+n] + A.fArray[ 4*N+n]*B.fArray[25*N+n] + A.fArray[ 5*N+n]*B.fArray[31*N+n];
      C.fArray[ 2*N+n] = A.fArray[ 0*N+n]*B.fArray[ 2*N+n] + A.fArray[ 1*N+n]*B.fArray[ 8*N+n] + A.fArray[ 2*N+n]*B.fArray[14*N+n] + A.fArray[ 3*N+n]*B.fArray[20*N+n] + A.fArray[ 4*N+n]*B.fArray[26*N+n] + A.fArray[ 5*N+n]*B.fArray[32*N+n];
      C.fArray[ 3*N+n] = A.fArray[ 0*N+n]*B.fArray[ 3*N+n] + A.fArray[ 1*N+n]*B.fArray[ 9*N+n] + A.fArray[ 2*N+n]*B.fArray[15*N+n] + A.fArray[ 3*N+n]*B.fArray[21*N+n] + A.fArray[ 4*N+n]*B.fArray[27*N+n] + A.fArray[ 5*N+n]*B.fArray[33*N+n];
      C.fArray[ 4*N+n] = A.fArray[ 0*N+n]*B.fArray[ 4*N+n] + A.fArray[ 1*N+n]*B.fArray[10*N+n] + A.fArray[ 2*N+n]*B.fArray[16*N+n] + A.fArray[ 3*N+n]*B.fArray[22*N+n] + A.fArray[ 4*N+n]*B.fArray[28*N+n] + A.fArray[ 5*N+n]*B.fArray[34*N+n];
      C.fArray[ 5*N+n] = A.fArray[ 0*N+n]*B.fArray[ 5*N+n] + A.fArray[ 1*N+n]*B.fArray[11*N+n] + A.fArray[ 2*N+n]*B.fArray[17*N+n] + A.fArray[ 3*N+n]*B.fArray[23*N+n] + A.fArray[ 4*N+n]*B.fArray[29*N+n] + A.fArray[ 5*N+n]*B.fArray[35*N+n];
      C.fArray[ 6*N+n] = A.fArray[ 6*N+n]*B.fArray[ 0*N+n] + A.fArray[ 7*N+n]*B.fArray[ 6*N+n] + A.fArray[ 8*N+n]*B.fArray[12*N+n] + A.fArray[ 9*N+n]*B.fArray[18*N+n] + A.fArray[10*N+n]*B.fArray[24*N+n] + A.fArray[11*N+n]*B.fArray[30*N+n];
      C.fArray[ 7*N+n] = A.fArray[ 6*N+n]*B.fArray[ 1*N+n] + A.fArray[ 7*N+n]*B.fArray[ 7*N+n] + A.fArray[ 8*N+n]*B.fArray[13*N+n] + A.fArray[ 9*N+n]*B.fArray[19*N+n] + A.fArray[10*N+n]*B.fArray[25*N+n] + A.fArray[11*N+n]*B.fArray[31*N+n];
      C.fArray[ 8*N+n] = A.fArray[ 6*N+n]*B.fArray[ 2*N+n] + A.fArray[ 7*N+n]*B.fArray[ 8*N+n] + A.fArray[ 8*N+n]*B.fArray[14*N+n] + A.fArray[ 9*N+n]*B.fArray[20*N+n] + A.fArray[10*N+n]*B.fArray[26*N+n] + A.fArray[11*N+n]*B.fArray[32*N+n];
      C.fArray[ 9*N+n] = A.fArray[ 6*N+n]*B.fArray[ 3*N+n] + A.fArray[ 7*N+n]*B.fArray[ 9*N+n] + A.fArray[ 8*N+n]*B.fArray[15*N+n] + A.fArray[ 9*N+n]*B.fArray[21*N+n] + A.fArray[10*N+n]*B.fArray[27*N+n] + A.fArray[11*N+n]*B.fArray[33*N+n];
      C.fArray[10*N+n] = A.fArray[ 6*N+n]*B.fArray[ 4*N+n] + A.fArray[ 7*N+n]*B.fArray[10*N+n] + A.fArray[ 8*N+n]*B.fArray[16*N+n] + A.fArray[ 9*N+n]*B.fArray[22*N+n] + A.fArray[10*N+n]*B.fArray[28*N+n] + A.fArray[11*N+n]*B.fArray[34*N+n];
      C.fArray[11*N+n] = A.fArray[ 6*N+n]*B.fArray[ 5*N+n] + A.fArray[ 7*N+n]*B.fArray[11*N+n] + A.fArray[ 8*N+n]*B.fArray[17*N+n] + A.fArray[ 9*N+n]*B.fArray[23*N+n] + A.fArray[10*N+n]*B.fArray[29*N+n] + A.fArray[11*N+n]*B.fArray[35*N+n];
      C.fArray[12*N+n] = A.fArray[12*N+n]*B.fArray[ 0*N+n] + A.fArray[13*N+n]*B.fArray[ 6*N+n] + A.fArray[14*N+n]*B.fArray[12*N+n] + A.fArray[15*N+n]*B.fArray[18*N+n] + A.fArray[16*N+n]*B.fArray[24*N+n] + A.fArray[17*N+n]*B.fArray[30*N+n];
      C.fArray[13*N+n] = A.fArray[12*N+n]*B.fArray[ 1*N+n] + A.fArray[13*N+n]*B.fArray[ 7*N+n] + A.fArray[14*N+n]*B.fArray[13*N+n] + A.fArray[15*N+n]*B.fArray[19*N+n] + A.fArray[16*N+n]*B.fArray[25*N+n] + A.fArray[17*N+n]*B.fArray[31*N+n];
      C.fArray[14*N+n] = A.fArray[12*N+n]*B.fArray[ 2*N+n] + A.fArray[13*N+n]*B.fArray[ 8*N+n] + A.fArray[14*N+n]*B.fArray[14*N+n] + A.fArray[15*N+n]*B.fArray[20*N+n] + A.fArray[16*N+n]*B.fArray[26*N+n] + A.fArray[17*N+n]*B.fArray[32*N+n];
      C.fArray[15*N+n] = A.fArray[12*N+n]*B.fArray[ 3*N+n] + A.fArray[13*N+n]*B.fArray[ 9*N+n] + A.fArray[14*N+n]*B.fArray[15*N+n] + A.fArray[15*N+n]*B.fArray[21*N+n] + A.fArray[16*N+n]*B.fArray[27*N+n] + A.fArray[17*N+n]*B.fArray[33*N+n];
      C.fArray[16*N+n] = A.fArray[12*N+n]*B.fArray[ 4*N+n] + A.fArray[13*N+n]*B.fArray[10*N+n] + A.fArray[14*N+n]*B.fArray[16*N+n] + A.fArray[15*N+n]*B.fArray[22*N+n] + A.fArray[16*N+n]*B.fArray[28*N+n] + A.fArray[17*N+n]*B.fArray[34*N+n];
      C.fArray[17*N+n] = A.fArray[12*N+n]*B.fArray[ 5*N+n] + A.fArray[13*N+n]*B.fArray[11*N+n] + A.fArray[14*N+n]*B.fArray[17*N+n] + A.fArray[15*N+n]*B.fArray[23*N+n] + A.fArray[16*N+n]*B.fArray[29*N+n] + A.fArray[17*N+n]*B.fArray[35*N+n];
      C.fArray[18*N+n] = A.fArray[18*N+n]*B.fArray[ 0*N+n] + A.fArray[19*N+n]*B.fArray[ 6*N+n] + A.fArray[20*N+n]*B.fArray[12*N+n] + A.fArray[21*N+n]*B.fArray[18*N+n] + A.fArray[22*N+n]*B.fArray[24*N+n] + A.fArray[23*N+n]*B.fArray[30*N+n];
      C.fArray[19*N+n] = A.fArray[18*N+n]*B.fArray[ 1*N+n] + A.fArray[19*N+n]*B.fArray[ 7*N+n] + A.fArray[20*N+n]*B.fArray[13*N+n] + A.fArray[21*N+n]*B.fArray[19*N+n] + A.fArray[22*N+n]*B.fArray[25*N+n] + A.fArray[23*N+n]*B.fArray[31*N+n];
      C.fArray[20*N+n] = A.fArray[18*N+n]*B.fArray[ 2*N+n] + A.fArray[19*N+n]*B.fArray[ 8*N+n] + A.fArray[20*N+n]*B.fArray[14*N+n] + A.fArray[21*N+n]*B.fArray[20*N+n] + A.fArray[22*N+n]*B.fArray[26*N+n] + A.fArray[23*N+n]*B.fArray[32*N+n];
      C.fArray[21*N+n] = A.fArray[18*N+n]*B.fArray[ 3*N+n] + A.fArray[19*N+n]*B.fArray[ 9*N+n] + A.fArray[20*N+n]*B.fArray[15*N+n] + A.fArray[21*N+n]*B.fArray[21*N+n] + A.fArray[22*N+n]*B.fArray[27*N+n] + A.fArray[23*N+n]*B.fArray[33*N+n];
      C.fArray[22*N+n] = A.fArray[18*N+n]*B.fArray[ 4*N+n] + A.fArray[19*N+n]*B.fArray[10*N+n] + A.fArray[20*N+n]*B.fArray[16*N+n] + A.fArray[21*N+n]*B.fArray[22*N+n] + A.fArray[22*N+n]*B.fArray[28*N+n] + A.fArray[23*N+n]*B.fArray[34*N+n];
      C.fArray[23*N+n] = A.fArray[18*N+n]*B.fArray[ 5*N+n] + A.fArray[19*N+n]*B.fArray[11*N+n] + A.fArray[20*N+n]*B.fArray[17*N+n] + A.fArray[21*N+n]*B.fArray[23*N+n] + A.fArray[22*N+n]*B.fArray[29*N+n] + A.fArray[23*N+n]*B.fArray[35*N+n];
      C.fArray[24*N+n] = A.fArray[24*N+n]*B.fArray[ 0*N+n] + A.fArray[25*N+n]*B.fArray[ 6*N+n] + A.fArray[26*N+n]*B.fArray[12*N+n] + A.fArray[27*N+n]*B.fArray[18*N+n] + A.fArray[28*N+n]*B.fArray[24*N+n] + A.fArray[29*N+n]*B.fArray[30*N+n];
      C.fArray[25*N+n] = A.fArray[24*N+n]*B.fArray[ 1*N+n] + A.fArray[25*N+n]*B.fArray[ 7*N+n] + A.fArray[26*N+n]*B.fArray[13*N+n] + A.fArray[27*N+n]*B.fArray[19*N+n] + A.fArray[28*N+n]*B.fArray[25*N+n] + A.fArray[29*N+n]*B.fArray[31*N+n];
      C.fArray[26*N+n] = A.fArray[24*N+n]*B.fArray[ 2*N+n] + A.fArray[25*N+n]*B.fArray[ 8*N+n] + A.fArray[26*N+n]*B.fArray[14*N+n] + A.fArray[27*N+n]*B.fArray[20*N+n] + A.fArray[28*N+n]*B.fArray[26*N+n] + A.fArray[29*N+n]*B.fArray[32*N+n];
      C.fArray[27*N+n] = A.fArray[24*N+n]*B.fArray[ 3*N+n] + A.fArray[25*N+n]*B.fArray[ 9*N+n] + A.fArray[26*N+n]*B.fArray[15*N+n] + A.fArray[27*N+n]*B.fArray[21*N+n] + A.fArray[28*N+n]*B.fArray[27*N+n] + A.fArray[29*N+n]*B.fArray[33*N+n];
      C.fArray[28*N+n] = A.fArray[24*N+n]*B.fArray[ 4*N+n] + A.fArray[25*N+n]*B.fArray[10*N+n] + A.fArray[26*N+n]*B.fArray[16*N+n] + A.fArray[27*N+n]*B.fArray[22*N+n] + A.fArray[28*N+n]*B.fArray[28*N+n] + A.fArray[29*N+n]*B.fArray[34*N+n];
      C.fArray[29*N+n] = A.fArray[24*N+n]*B.fArray[ 5*N+n] + A.fArray[25*N+n]*B.fArray[11*N+n] + A.fArray[26*N+n]*B.fArray[17*N+n] + A.fArray[27*N+n]*B.fArray[23*N+n] + A.fArray[28*N+n]*B.fArray[29*N+n] + A.fArray[29*N+n]*B.fArray[35*N+n];
      C.fArray[30*N+n] = A.fArray[30*N+n]*B.fArray[ 0*N+n] + A.fArray[31*N+n]*B.fArray[ 6*N+n] + A.fArray[32*N+n]*B.fArray[12*N+n] + A.fArray[33*N+n]*B.fArray[18*N+n] + A.fArray[34*N+n]*B.fArray[24*N+n] + A.fArray[35*N+n]*B.fArray[30*N+n];
      C.fArray[31*N+n] = A.fArray[30*N+n]*B.fArray[ 1*N+n] + A.fArray[31*N+n]*B.fArray[ 7*N+n] + A.fArray[32*N+n]*B.fArray[13*N+n] + A.fArray[33*N+n]*B.fArray[19*N+n] + A.fArray[34*N+n]*B.fArray[25*N+n] + A.fArray[35*N+n]*B.fArray[31*N+n];
      C.fArray[32*N+n] = A.fArray[30*N+n]*B.fArray[ 2*N+n] + A.fArray[31*N+n]*B.fArray[ 8*N+n] + A.fArray[32*N+n]*B.fArray[14*N+n] + A.fArray[33*N+n]*B.fArray[20*N+n] + A.fArray[34*N+n]*B.fArray[26*N+n] + A.fArray[35*N+n]*B.fArray[32*N+n];
      C.fArray[33*N+n] = A.fArray[30*N+n]*B.fArray[ 3*N+n] + A.fArray[31*N+n]*B.fArray[ 9*N+n] + A.fArray[32*N+n]*B.fArray[15*N+n] + A.fArray[33*N+n]*B.fArray[21*N+n] + A.fArray[34*N+n]*B.fArray[27*N+n] + A.fArray[35*N+n]*B.fArray[33*N+n];
      C.fArray[34*N+n] = A.fArray[30*N+n]*B.fArray[ 4*N+n] + A.fArray[31*N+n]*B.fArray[10*N+n] + A.fArray[32*N+n]*B.fArray[16*N+n] + A.fArray[33*N+n]*B.fArray[22*N+n] + A.fArray[34*N+n]*B.fArray[28*N+n] + A.fArray[35*N+n]*B.fArray[34*N+n];
      C.fArray[35*N+n] = A.fArray[30*N+n]*B.fArray[ 5*N+n] + A.fArray[31*N+n]*B.fArray[11*N+n] + A.fArray[32*N+n]*B.fArray[17*N+n] + A.fArray[33*N+n]*B.fArray[23*N+n] + A.fArray[34*N+n]*B.fArray[29*N+n] + A.fArray[35*N+n]*B.fArray[35*N+n];
   }
}


//==============================================================================
// Cramer inversion
//==============================================================================

template<typename T, idx_t D, idx_t N>
struct CramerInverter
{
   static void Invert(MPlex<T, D, D, N>& A, double *determ=0)
   {
      throw std::runtime_error("general cramer inversion not supported");
   }
};


template<typename T, idx_t N>
struct CramerInverter<T, 2, N>
{
   static void Invert(MPlex<T, 2, 2, N>& A, double *determ=0)
   {
      typedef T TT;

      T *M = A.fArray;

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         const TT det = M[n] * M[3*N + n] - M[2*N + n] * M[N + n];

         //if (determ)
         //determ[n] = s;

         //if (det == 0)
         {
            const TT s   = TT(1) / det;
            const TT tmp = s * M[3*N + n];
            M[N + n]   *= -s;
            M[2*N + n] *= -s;
            M[3*N + n]  = s * M[n];
            M[n] = tmp;
         }
      }
   }
};

template<typename T, idx_t N>
struct CramerInverter<T, 3, N>
{
   static void Invert(MPlex<T, 3, 3, N>& A, double *determ=0)
   {
      typedef T TT;

      T *M = A.fArray;

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         const TT c00 = M[4*N+n] * M[8*N+n] - M[5*N+n] * M[7*N+n];
         const TT c01 = M[5*N+n] * M[6*N+n] - M[3*N+n] * M[8*N+n];
         const TT c02 = M[3*N+n] * M[7*N+n] - M[4*N+n] * M[6*N+n];
         const TT c10 = M[7*N+n] * M[2*N+n] - M[8*N+n] * M[1*N+n];
         const TT c11 = M[8*N+n] * M[0*N+n] - M[6*N+n] * M[2*N+n];
         const TT c12 = M[6*N+n] * M[1*N+n] - M[7*N+n] * M[0*N+n];
         const TT c20 = M[1*N+n] * M[5*N+n] - M[2*N+n] * M[4*N+n];
         const TT c21 = M[2*N+n] * M[3*N+n] - M[0*N+n] * M[5*N+n];
         const TT c22 = M[0*N+n] * M[4*N+n] - M[1*N+n] * M[3*N+n];

         const TT det = M[0*N+n] * c00 + M[1*N+n] * c01 + M[2*N+n] * c02;

         //if (determ)
         //  *determ[n] = det;

         const TT s = TT(1) / det;

         M[0*N+n] = s*c00;
         M[1*N+n] = s*c10;
         M[2*N+n] = s*c20;
         M[3*N+n] = s*c01;
         M[4*N+n] = s*c11;
         M[5*N+n] = s*c21;
         M[6*N+n] = s*c02;
         M[7*N+n] = s*c12;
         M[8*N+n] = s*c22;
      }
   }
};

template<typename T, idx_t D, idx_t N>
void InvertCramer(MPlex<T, D, D, N>& A, double *determ=0)
{
   // We don't do general Inverts.

   CramerInverter<T, D, N>::Invert(A, determ);
}


//==============================================================================
// Cholesky inversion
//==============================================================================

template<typename T, idx_t D, idx_t N>
struct CholeskyInverter
{
   static void Invert(MPlex<T, D, D, N>& A)
   {
      throw std::runtime_error("general cholesky inversion not supported");
   }
};

template<typename T, idx_t N>
struct CholeskyInverter<T, 3, N>
{
   // Remember, this only works on symmetric matrices!
   // Optimized version for positive definite matrices, no checks.
   // Also, use as little locals as possible.
   // This gives: host  x 5.8 (instead of 4.7x)
   //             mic   x17.7 (instead of 8.5x))
   static void Invert(MPlex<T, 3, 3, N>& A)
   {
      typedef T TT;

      T *M = A.fArray;

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         TT l0 = std::sqrt(T(1) / M[0*N+n]);
         TT l1 = M[3*N+n] * l0;
         TT l2 = M[4*N+n] - l1 * l1;
         l2 = std::sqrt(T(1) / l2);
         TT l3 = M[6*N+n] * l0;
         TT l4 = (M[7*N+n] - l1 * l3) * l2;
         TT l5 = M[8*N+n] - (l3 * l3 + l4 * l4);
         l5 = std::sqrt(T(1) / l5);

         // decomposition done

         l3 = (l1 * l4 * l2 - l3) * l0 * l5;
         l1 = -l1 * l0 * l2;
         l4 = -l4 * l2 * l5;

         M[0*N+n] = l3*l3 + l1*l1 + l0*l0;
         M[1*N+n] = M[3*N+n] = l3*l4 + l1*l2;
         M[4*N+n] = l4*l4 + l2*l2;
         M[2*N+n] = M[6*N+n] = l3*l5;
         M[5*N+n] = M[7*N+n] = l4*l5;
         M[8*N+n] = l5*l5;

         // m(2,x) are all zero if anything went wrong at l5.
         // all zero, if anything went wrong already for l0 or l2.
      }
   }
};

template<typename T, idx_t D, idx_t N>
void InvertCholesky(MPlex<T, D, D, N>& A)
{
   CholeskyInverter<T, D, N>::Invert(A);
}

}

#endif
