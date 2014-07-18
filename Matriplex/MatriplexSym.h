#ifndef MatriplexSym_H
#define MatriplexSym_H

#include "MatriplexCommon.h"
#include "Matriplex.h"

namespace Matriplex
{

  const idx_t gSymOffsets[7][36] =
  {
    {}, {},
    { 0, 1, 1, 2 },
    { 0, 1, 3, 1, 2, 4, 3, 4, 5 }, // 3
    {}, {},
    { 0, 1, 3, 6, 10, 15, 1, 2, 4, 7, 11, 16, 3, 4, 5, 8, 12, 17, 6, 7, 8, 9, 13, 18, 10, 11, 12, 13, 14, 19, 15, 16, 17, 18, 19, 20 }
  };

//------------------------------------------------------------------------------

template<typename T, idx_t D, idx_t N>
class MatriplexSym
{
public:
   typedef T value_type;

   enum
   {
      /// no. of matrix rows
      kRows = D,
      /// no. of matrix columns
      kCols = D,
      /// no of elements: lower triangle
      kSize = (D + 1) * D / 2,
      /// size of the whole matriplex
      kTotSize = N * kSize
   };

   T fArray[kTotSize] __attribute__((aligned(64)));


   MatriplexSym()    {}
   MatriplexSym(T v) { SetVal(v); }

   void SetVal(T v)
   {
      for (idx_t i = 0; i < kTotSize; ++i)
      {
         fArray[i] = v;
      }
   }

   idx_t * Offsets()    { return gSymOffsets[D];    }
   idx_t   Off(idx_t i) { return gSymOffsets[D][i]; }

   // This is krappe
   T& At(idx_t n, idx_t i, idx_t j) { return fArray[Off(i * D + j) * N + n]; }

   T& operator()(idx_t n, idx_t i, idx_t j) { return At(n, i, j); }

   void Assign(idx_t n, T *arr)
   {
#pragma simd
      for (idx_t i = n; i < kTotSize; i += N)
      {
         fArray[i] = *(arr++);
      }
   }
};


template<typename T, idx_t D, idx_t N> using MPlexSym = MatriplexSym<T, D, N>;


//==============================================================================
// Multiplications
//==============================================================================

template<typename T, idx_t D, idx_t N>
struct SymMultiplyCls
{
   static void Multiply(const MatriplexSym<T, D, N>& A,
                        const MatriplexSym<T, D, N>& B,
                        Matriplex<T, D, D, N>& C)
   {
      throw std::runtime_error("general symmetric multiplication not supported");
   }
};


template<typename T, idx_t N>
struct SymMultiplyCls<T, 3, N>
{
   static void Multiply(const MatriplexSym<T, 3, N>& A,
                        const MatriplexSym<T, 3, N>& B,
                        Matriplex<T, 3, 3, N>& C)
{
   const T *a = A.fArray, *b = B.fArray;
         T *c = C.fArray;

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      c[ 0*N+n] = a[ 0*N+n]*b[ 0*N+n] + a[ 1*N+n]*b[ 1*N+n] + a[ 3*N+n]*b[ 3*N+n];
      c[ 1*N+n] = a[ 0*N+n]*b[ 1*N+n] + a[ 1*N+n]*b[ 2*N+n] + a[ 3*N+n]*b[ 4*N+n];
      c[ 2*N+n] = a[ 0*N+n]*b[ 3*N+n] + a[ 1*N+n]*b[ 4*N+n] + a[ 3*N+n]*b[ 5*N+n];
      c[ 3*N+n] = a[ 1*N+n]*b[ 0*N+n] + a[ 2*N+n]*b[ 1*N+n] + a[ 4*N+n]*b[ 3*N+n];
      c[ 4*N+n] = a[ 1*N+n]*b[ 1*N+n] + a[ 2*N+n]*b[ 2*N+n] + a[ 4*N+n]*b[ 4*N+n];
      c[ 5*N+n] = a[ 1*N+n]*b[ 3*N+n] + a[ 2*N+n]*b[ 4*N+n] + a[ 4*N+n]*b[ 5*N+n];
      c[ 6*N+n] = a[ 3*N+n]*b[ 0*N+n] + a[ 4*N+n]*b[ 1*N+n] + a[ 5*N+n]*b[ 3*N+n];
      c[ 7*N+n] = a[ 3*N+n]*b[ 1*N+n] + a[ 4*N+n]*b[ 2*N+n] + a[ 5*N+n]*b[ 4*N+n];
      c[ 8*N+n] = a[ 3*N+n]*b[ 3*N+n] + a[ 4*N+n]*b[ 4*N+n] + a[ 5*N+n]*b[ 5*N+n];
   }
}
};

template<typename T, idx_t N>
struct SymMultiplyCls<T, 6, N>
{
   static void Multiply(const MatriplexSym<float, 6, N>& A,
                        const MatriplexSym<float, 6, N>& B,
                        Matriplex<float, 6, 6, N>& C)
{
   const T *a = A.fArray, *b = B.fArray;
         T *c = C.fArray;

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      c[ 0*N+n] = a[ 0*N+n]*b[ 0*N+n] + a[ 1*N+n]*b[ 1*N+n] + a[ 3*N+n]*b[ 3*N+n] + a[ 6*N+n]*b[ 6*N+n] + a[10*N+n]*b[10*N+n] + a[15*N+n]*b[15*N+n];
      c[ 1*N+n] = a[ 0*N+n]*b[ 1*N+n] + a[ 1*N+n]*b[ 2*N+n] + a[ 3*N+n]*b[ 4*N+n] + a[ 6*N+n]*b[ 7*N+n] + a[10*N+n]*b[11*N+n] + a[15*N+n]*b[16*N+n];
      c[ 2*N+n] = a[ 0*N+n]*b[ 3*N+n] + a[ 1*N+n]*b[ 4*N+n] + a[ 3*N+n]*b[ 5*N+n] + a[ 6*N+n]*b[ 8*N+n] + a[10*N+n]*b[12*N+n] + a[15*N+n]*b[17*N+n];
      c[ 3*N+n] = a[ 0*N+n]*b[ 6*N+n] + a[ 1*N+n]*b[ 7*N+n] + a[ 3*N+n]*b[ 8*N+n] + a[ 6*N+n]*b[ 9*N+n] + a[10*N+n]*b[13*N+n] + a[15*N+n]*b[18*N+n];
      c[ 4*N+n] = a[ 0*N+n]*b[10*N+n] + a[ 1*N+n]*b[11*N+n] + a[ 3*N+n]*b[12*N+n] + a[ 6*N+n]*b[13*N+n] + a[10*N+n]*b[14*N+n] + a[15*N+n]*b[19*N+n];
      c[ 5*N+n] = a[ 0*N+n]*b[15*N+n] + a[ 1*N+n]*b[16*N+n] + a[ 3*N+n]*b[17*N+n] + a[ 6*N+n]*b[18*N+n] + a[10*N+n]*b[19*N+n] + a[15*N+n]*b[20*N+n];
      c[ 6*N+n] = a[ 1*N+n]*b[ 0*N+n] + a[ 2*N+n]*b[ 1*N+n] + a[ 4*N+n]*b[ 3*N+n] + a[ 7*N+n]*b[ 6*N+n] + a[11*N+n]*b[10*N+n] + a[16*N+n]*b[15*N+n];
      c[ 7*N+n] = a[ 1*N+n]*b[ 1*N+n] + a[ 2*N+n]*b[ 2*N+n] + a[ 4*N+n]*b[ 4*N+n] + a[ 7*N+n]*b[ 7*N+n] + a[11*N+n]*b[11*N+n] + a[16*N+n]*b[16*N+n];
      c[ 8*N+n] = a[ 1*N+n]*b[ 3*N+n] + a[ 2*N+n]*b[ 4*N+n] + a[ 4*N+n]*b[ 5*N+n] + a[ 7*N+n]*b[ 8*N+n] + a[11*N+n]*b[12*N+n] + a[16*N+n]*b[17*N+n];
      c[ 9*N+n] = a[ 1*N+n]*b[ 6*N+n] + a[ 2*N+n]*b[ 7*N+n] + a[ 4*N+n]*b[ 8*N+n] + a[ 7*N+n]*b[ 9*N+n] + a[11*N+n]*b[13*N+n] + a[16*N+n]*b[18*N+n];
      c[10*N+n] = a[ 1*N+n]*b[10*N+n] + a[ 2*N+n]*b[11*N+n] + a[ 4*N+n]*b[12*N+n] + a[ 7*N+n]*b[13*N+n] + a[11*N+n]*b[14*N+n] + a[16*N+n]*b[19*N+n];
      c[11*N+n] = a[ 1*N+n]*b[15*N+n] + a[ 2*N+n]*b[16*N+n] + a[ 4*N+n]*b[17*N+n] + a[ 7*N+n]*b[18*N+n] + a[11*N+n]*b[19*N+n] + a[16*N+n]*b[20*N+n];
      c[12*N+n] = a[ 3*N+n]*b[ 0*N+n] + a[ 4*N+n]*b[ 1*N+n] + a[ 5*N+n]*b[ 3*N+n] + a[ 8*N+n]*b[ 6*N+n] + a[12*N+n]*b[10*N+n] + a[17*N+n]*b[15*N+n];
      c[13*N+n] = a[ 3*N+n]*b[ 1*N+n] + a[ 4*N+n]*b[ 2*N+n] + a[ 5*N+n]*b[ 4*N+n] + a[ 8*N+n]*b[ 7*N+n] + a[12*N+n]*b[11*N+n] + a[17*N+n]*b[16*N+n];
      c[14*N+n] = a[ 3*N+n]*b[ 3*N+n] + a[ 4*N+n]*b[ 4*N+n] + a[ 5*N+n]*b[ 5*N+n] + a[ 8*N+n]*b[ 8*N+n] + a[12*N+n]*b[12*N+n] + a[17*N+n]*b[17*N+n];
      c[15*N+n] = a[ 3*N+n]*b[ 6*N+n] + a[ 4*N+n]*b[ 7*N+n] + a[ 5*N+n]*b[ 8*N+n] + a[ 8*N+n]*b[ 9*N+n] + a[12*N+n]*b[13*N+n] + a[17*N+n]*b[18*N+n];
      c[16*N+n] = a[ 3*N+n]*b[10*N+n] + a[ 4*N+n]*b[11*N+n] + a[ 5*N+n]*b[12*N+n] + a[ 8*N+n]*b[13*N+n] + a[12*N+n]*b[14*N+n] + a[17*N+n]*b[19*N+n];
      c[17*N+n] = a[ 3*N+n]*b[15*N+n] + a[ 4*N+n]*b[16*N+n] + a[ 5*N+n]*b[17*N+n] + a[ 8*N+n]*b[18*N+n] + a[12*N+n]*b[19*N+n] + a[17*N+n]*b[20*N+n];
      c[18*N+n] = a[ 6*N+n]*b[ 0*N+n] + a[ 7*N+n]*b[ 1*N+n] + a[ 8*N+n]*b[ 3*N+n] + a[ 9*N+n]*b[ 6*N+n] + a[13*N+n]*b[10*N+n] + a[18*N+n]*b[15*N+n];
      c[19*N+n] = a[ 6*N+n]*b[ 1*N+n] + a[ 7*N+n]*b[ 2*N+n] + a[ 8*N+n]*b[ 4*N+n] + a[ 9*N+n]*b[ 7*N+n] + a[13*N+n]*b[11*N+n] + a[18*N+n]*b[16*N+n];
      c[20*N+n] = a[ 6*N+n]*b[ 3*N+n] + a[ 7*N+n]*b[ 4*N+n] + a[ 8*N+n]*b[ 5*N+n] + a[ 9*N+n]*b[ 8*N+n] + a[13*N+n]*b[12*N+n] + a[18*N+n]*b[17*N+n];
      c[21*N+n] = a[ 6*N+n]*b[ 6*N+n] + a[ 7*N+n]*b[ 7*N+n] + a[ 8*N+n]*b[ 8*N+n] + a[ 9*N+n]*b[ 9*N+n] + a[13*N+n]*b[13*N+n] + a[18*N+n]*b[18*N+n];
      c[22*N+n] = a[ 6*N+n]*b[10*N+n] + a[ 7*N+n]*b[11*N+n] + a[ 8*N+n]*b[12*N+n] + a[ 9*N+n]*b[13*N+n] + a[13*N+n]*b[14*N+n] + a[18*N+n]*b[19*N+n];
      c[23*N+n] = a[ 6*N+n]*b[15*N+n] + a[ 7*N+n]*b[16*N+n] + a[ 8*N+n]*b[17*N+n] + a[ 9*N+n]*b[18*N+n] + a[13*N+n]*b[19*N+n] + a[18*N+n]*b[20*N+n];
      c[24*N+n] = a[10*N+n]*b[ 0*N+n] + a[11*N+n]*b[ 1*N+n] + a[12*N+n]*b[ 3*N+n] + a[13*N+n]*b[ 6*N+n] + a[14*N+n]*b[10*N+n] + a[19*N+n]*b[15*N+n];
      c[25*N+n] = a[10*N+n]*b[ 1*N+n] + a[11*N+n]*b[ 2*N+n] + a[12*N+n]*b[ 4*N+n] + a[13*N+n]*b[ 7*N+n] + a[14*N+n]*b[11*N+n] + a[19*N+n]*b[16*N+n];
      c[26*N+n] = a[10*N+n]*b[ 3*N+n] + a[11*N+n]*b[ 4*N+n] + a[12*N+n]*b[ 5*N+n] + a[13*N+n]*b[ 8*N+n] + a[14*N+n]*b[12*N+n] + a[19*N+n]*b[17*N+n];
      c[27*N+n] = a[10*N+n]*b[ 6*N+n] + a[11*N+n]*b[ 7*N+n] + a[12*N+n]*b[ 8*N+n] + a[13*N+n]*b[ 9*N+n] + a[14*N+n]*b[13*N+n] + a[19*N+n]*b[18*N+n];
      c[28*N+n] = a[10*N+n]*b[10*N+n] + a[11*N+n]*b[11*N+n] + a[12*N+n]*b[12*N+n] + a[13*N+n]*b[13*N+n] + a[14*N+n]*b[14*N+n] + a[19*N+n]*b[19*N+n];
      c[29*N+n] = a[10*N+n]*b[15*N+n] + a[11*N+n]*b[16*N+n] + a[12*N+n]*b[17*N+n] + a[13*N+n]*b[18*N+n] + a[14*N+n]*b[19*N+n] + a[19*N+n]*b[20*N+n];
      c[30*N+n] = a[15*N+n]*b[ 0*N+n] + a[16*N+n]*b[ 1*N+n] + a[17*N+n]*b[ 3*N+n] + a[18*N+n]*b[ 6*N+n] + a[19*N+n]*b[10*N+n] + a[20*N+n]*b[15*N+n];
      c[31*N+n] = a[15*N+n]*b[ 1*N+n] + a[16*N+n]*b[ 2*N+n] + a[17*N+n]*b[ 4*N+n] + a[18*N+n]*b[ 7*N+n] + a[19*N+n]*b[11*N+n] + a[20*N+n]*b[16*N+n];
      c[32*N+n] = a[15*N+n]*b[ 3*N+n] + a[16*N+n]*b[ 4*N+n] + a[17*N+n]*b[ 5*N+n] + a[18*N+n]*b[ 8*N+n] + a[19*N+n]*b[12*N+n] + a[20*N+n]*b[17*N+n];
      c[33*N+n] = a[15*N+n]*b[ 6*N+n] + a[16*N+n]*b[ 7*N+n] + a[17*N+n]*b[ 8*N+n] + a[18*N+n]*b[ 9*N+n] + a[19*N+n]*b[13*N+n] + a[20*N+n]*b[18*N+n];
      c[34*N+n] = a[15*N+n]*b[10*N+n] + a[16*N+n]*b[11*N+n] + a[17*N+n]*b[12*N+n] + a[18*N+n]*b[13*N+n] + a[19*N+n]*b[14*N+n] + a[20*N+n]*b[19*N+n];
      c[35*N+n] = a[15*N+n]*b[15*N+n] + a[16*N+n]*b[16*N+n] + a[17*N+n]*b[17*N+n] + a[18*N+n]*b[18*N+n] + a[19*N+n]*b[19*N+n] + a[20*N+n]*b[20*N+n];
   }
}
};

template<typename T, idx_t D, idx_t N>
void Multiply(const MatriplexSym<T, D, N>& A,
              const MatriplexSym<T, D, N>& B,
                    Matriplex<T, D, D, N>& C)
{
   // printf("Multipl %d %d\n", D, N);

   SymMultiplyCls<T, D, N>::Multiply(A, B, C);
}


//==============================================================================
// Cramer inversion
//==============================================================================

template<typename T, idx_t D, idx_t N>
struct CramerInverterSym
{
   static void Invert(MatriplexSym<T, D, N>& A)
   {
     throw std::runtime_error("general cramer inversion not supported");
   }
};

template<typename T, idx_t N>
struct CramerInverterSym<T, 2, N>
{
   static void Invert(MatriplexSym<T, 2, N>& A)
   {
      typedef T TT;

      T *a = A.fArray;

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         const TT det = a[0*N+n] * a[2*N+n] - a[1*N+n] * a[1*N+n];

         const TT s   = TT(1) / det;
         const TT tmp = s * a[2*N+n];
         a[1*N+n] *= -s;
         a[2*N+n]  = s * a[0*N+n];
         a[0*N+n]  = tmp;
      }
   }
};

template<typename T, idx_t N>
struct CramerInverterSym<T, 3, N>
{
   static void Invert(MatriplexSym<T, 3, N>& A)
   {
      typedef T TT;

      T *a = A.fArray;

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         const TT c00 = a[2*N+n] * a[5*N+n] - a[4*N+n] * a[4*N+n];
         const TT c01 = a[4*N+n] * a[3*N+n] - a[1*N+n] * a[5*N+n];
         const TT c02 = a[1*N+n] * a[4*N+n] - a[2*N+n] * a[3*N+n];
         const TT c11 = a[5*N+n] * a[0*N+n] - a[3*N+n] * a[3*N+n];
         const TT c12 = a[3*N+n] * a[1*N+n] - a[4*N+n] * a[0*N+n];
         const TT c22 = a[0*N+n] * a[2*N+n] - a[1*N+n] * a[1*N+n];

         const TT det = a[0*N+n] * c00 + a[1*N+n] * c01 + a[3*N+n] * c02;

         const TT s = TT(1) / det;

         a[0*N+n] = s*c00;
         a[1*N+n] = s*c01;
         a[2*N+n] = s*c11;
         a[3*N+n] = s*c02;
         a[4*N+n] = s*c12;
         a[5*N+n] = s*c22;
      }
   }
};

template<typename T, idx_t D, idx_t N>
void InvertCramerSym(MatriplexSym<T, D, N>& A)
{
   CramerInverterSym<T, D, N>::Invert(A);
}


//==============================================================================
// Cholesky inversion
//==============================================================================

template<typename T, idx_t D, idx_t N>
struct CholeskyInverterSym
{
   static void Invert(MatriplexSym<T, D, N>& A)
   {
     throw std::runtime_error("general cholesky inversion not supported");
   }
};

template<typename T, idx_t N>
struct CholeskyInverterSym<T, 3, N>
{
   static void Invert(MatriplexSym<T, 3, N>& A)
   {
      typedef T TT;

      T *a = A.fArray;

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         TT l0 = std::sqrt(T(1) / a[0*N+n]);
         TT l1 = a[1*N+n] * l0;
         TT l2 = a[2*N+n] - l1 * l1;
         l2 = std::sqrt(T(1) / l2);
         TT l3 = a[3*N+n] * l0;
         TT l4 = (a[4*N+n] - l1 * l3) * l2;
         TT l5 = a[5*N+n] - (l3 * l3 + l4 * l4);
         l5 = std::sqrt(T(1) / l5);

         // decomposition done

         l3 = (l1 * l4 * l2 - l3) * l0 * l5;
         l1 = -l1 * l0 * l2;
         l4 = -l4 * l2 * l5;

         a[0*N+n] = l3*l3 + l1*l1 + l0*l0;
         a[1*N+n] = l3*l4 + l1*l2;
         a[2*N+n] = l4*l4 + l2*l2;
         a[3*N+n] = l3*l5;
         a[4*N+n] = l4*l5;
         a[5*N+n] = l5*l5;

         // m(2,x) are all zero if anything went wrong at l5.
         // all zero, if anything went wrong already for l0 or l2.
      }
   }
};

template<typename T, idx_t D, idx_t N>
void InvertCholeskySym(MatriplexSym<T, D, N>& A)
{
   CholeskyInverterSym<T, D, N>::Invert(A);
}

}

#endif
