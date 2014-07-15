#ifndef Matriplex_H
#define Matriplex_H

#include "MatriplexCommon.h"

namespace Matriplex
{

template<typename T, idx_t D1, idx_t D2, idx_t N>
class Matriplex
{
public:

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


template<typename T, idx_t D1, idx_t D2, idx_t D3, idx_t N>
void Multiply(const Matriplex<T, D1, D2, N>& A,
              const Matriplex<T, D2, D3, N>& B,
              Matriplex<T, D1, D3, N>& C)
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
void MultiplyUnrolled(const Matriplex<T, D1, D2, N>& A,
                      const Matriplex<T, D2, D3, N>& B,
                      Matriplex<T, D1, D3, N>& C)
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
   CramerInverter(Matriplex<T, D, D, N>& C, double *determ=0)
   {
      throw std::runtime_error("general cramer inversion not supported");
   }
};


template<typename T, idx_t N>
struct CramerInverter<T, 2, N>
{
   void operator()(Matriplex<T, 2, 2, N>& C, double *determ=0)
   {
      typedef T TT;

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         T *pM = & C.fArray[0];

         const TT det = pM[n] * pM[3*N + n] - pM[2*N + n] * pM[N + n];

         //if (determ)
         //determ[n] = s;

         //if (det == 0)
         {
            const TT s   = TT(1) / det;
            const TT tmp = s * pM[3*N + n];
            pM[N + n]   *= -s;
            pM[2*N + n] *= -s;
            pM[3*N + n]  = s * pM[n];
            pM[n] = tmp;
         }
      }
   }
};

template<typename T, idx_t N>
struct CramerInverter<T, 3, N>
{
   void operator()(Matriplex<T, 3, 3, N>& C, double *determ=0)
   {
      typedef T TT;

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         T *pM = & C.fArray[n];

         const TT c00 = pM[4*N] * pM[8*N] - pM[5*N] * pM[7*N];
         const TT c01 = pM[5*N] * pM[6*N] - pM[3*N] * pM[8*N];
         const TT c02 = pM[3*N] * pM[7*N] - pM[4*N] * pM[6*N];
         const TT c10 = pM[7*N] * pM[2*N] - pM[8*N] * pM[1*N];
         const TT c11 = pM[8*N] * pM[0*N] - pM[6*N] * pM[2*N];
         const TT c12 = pM[6*N] * pM[1*N] - pM[7*N] * pM[0*N];
         const TT c20 = pM[1*N] * pM[5*N] - pM[2*N] * pM[4*N];
         const TT c21 = pM[2*N] * pM[3*N] - pM[0*N] * pM[5*N];
         const TT c22 = pM[0*N] * pM[4*N] - pM[1*N] * pM[3*N];

         const TT det = pM[0*N] * c00 + pM[1*N] * c01 + pM[2*N] * c02;

         //if (determ)
         //  *determ[n] = det;

         const TT s = TT(1) / det;

         pM[0*N] = s*c00;
         pM[1*N] = s*c10;
         pM[2*N] = s*c20;
         pM[3*N] = s*c01;
         pM[4*N] = s*c11;
         pM[5*N] = s*c21;
         pM[6*N] = s*c02;
         pM[7*N] = s*c12;
         pM[8*N] = s*c22;
      }
   }
};

template<typename T, idx_t D, idx_t N>
void InvertCramer(Matriplex<T, D, D, N>& C, double *determ=0)
{
   // We don't do general Inverts.

   CramerInverter<T, D, N> ci;
   ci(C, determ);
}


//==============================================================================
// Cholesky inversion
//==============================================================================

template<typename T, idx_t D, idx_t N>
struct CholInverter
{
   CholInverter(Matriplex<T, D, D, N>& m)
   {
      throw std::runtime_error("general cholesky inversion not supported");
   }
};


template<typename T, idx_t N>
struct CholInverter<T, 3, N>
{
   enum
   {
      kSize    = 3 * (3 + 1) / 2,
      kTotSize = N * kSize
   };

   //T  fL[kTotSize];

   //T& l(idx_t i, idx_t n) { return fL[i * kSize + n]; }

   /*
   // The "slow" version that does >=0 checks.
   CholInverter(Matriplex<T, 3, 3, N>& m)
   {
   #pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
   T l0 = (m(0,0,n) > T(0)) ? std::sqrt(T(1) / m(0,0,n)) : 0;
   T l1 = m(1,0,n) * l0;
   T l2 = m(1,1,n) - l1 * l1;
   l2 = (l2 > T(0)) ? std::sqrt(T(1) / l2) : 0;
   T l3 = m(2,0,n) * l0;
   T l4 = (m(2,1,n) - l1 * l3) * l2;
   T l5 = m(2,2,n) - (l3 * l3 + l4 * l4);
   l5 = (l5 > T(0)) ? std::sqrt(T(1) / l5) : 0;

   // decomposition done

   const T li21 = -l1 * l0 * l2;
   const T li32 = -l4 * l2 * l5;
   const T li31 = (l1 * l4 * l2 - l3) * l0 * l5;

   m(0,0,n) = li31*li31 + li21*li21 + l0*l0;
   m(1,0,n) = m(0,1,n) = li31*li32 + li21*l2;
   m(1,1,n) = li32*li32 + l2*l2;
   m(2,0,n) = m(0,2,n) = li31*l5;
   m(2,1,n) = m(1,2,n) = li32*l5;
   m(2,2,n) = l5*l5;

   // m(2,x) are all zero if anything went wrong at l5.
   // all zero, if anything went wrong already for l0 or l2.
   }
   }
   */

   /*
   // Optimized version for positive definite matrices, no checks.
   // Also, use as little locals as possible.
   // Fill only part of output matrix --> need MatriplexSym !!!
   // This gives: host  x 5.8 (instead of 4.7x)
   //             mic   x17.7 (instead of 8.5x))
   */
   CholInverter(Matriplex<T, 3, 3, N>& m)
   {
#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         T l0 = std::sqrt(T(1) / m(0,0,n));
         T l1 = m(1,0,n) * l0;
         T l2 = m(1,1,n) - l1 * l1;
         l2 = std::sqrt(T(1) / l2);
         T l3 = m(2,0,n) * l0;
         T l4 = (m(2,1,n) - l1 * l3) * l2;
         T l5 = m(2,2,n) - (l3 * l3 + l4 * l4);
         l5 = std::sqrt(T(1) / l5);

         // decomposition done

         l3 = (l1 * l4 * l2 - l3) * l0 * l5;
         l1 = -l1 * l0 * l2;
         l4 = -l4 * l2 * l5;

         m(0,0,n) = l3*l3 + l1*l1 + l0*l0;
         m(1,0,n) = l3*l4 + l1*l2;
         m(1,1,n) = l4*l4 + l2*l2;
         m(2,0,n) = l3*l5;
         m(2,1,n) = l4*l5;
         m(2,2,n) = l5*l5;

         // m(2,x) are all zero if anything went wrong at l5.
         // all zero, if anything went wrong already for l0 or l2.
      }
   }
};

template<typename T, idx_t D, idx_t N>
void InvertChol(Matriplex<T, D, D, N>& m)
{
   // We don't do general Inverts.

   CholInverter<T, D, N> ci(m);
   //ci(C, determ);
}

}

#endif
