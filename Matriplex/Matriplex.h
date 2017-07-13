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

   Matriplex()    {}
   Matriplex(T v) { SetVal(v); }

   idx_t PlexSize() const { return N; }

   void SetVal(T v)
   {
      for (idx_t i = 0; i < kTotSize; ++i)
      {
         fArray[i] = v;
      }
   }

   T  operator[](idx_t xx) const { return fArray[xx]; }
   T& operator[](idx_t xx)       { return fArray[xx]; }

   const T& ConstAt(idx_t n, idx_t i, idx_t j) const { return fArray[(i * D2 + j) * N + n]; }

   T& At(idx_t n, idx_t i, idx_t j) { return fArray[(i * D2 + j) * N + n]; }

   T& operator()(idx_t n, idx_t i, idx_t j) { return fArray[(i * D2 + j) * N + n]; }
   const T& operator()(idx_t n, idx_t i, idx_t j) const { return fArray[(i * D2 + j) * N + n]; }

   Matriplex& operator=(const Matriplex& m)
   {
      memcpy(fArray, m.fArray, sizeof(T) * kTotSize); return *this;
   }

   void CopyIn(idx_t n, const T *arr)
   {
      for (idx_t i = n; i < kTotSize; i += N)
      {
         fArray[i] = *(arr++);
      }
   }

#if defined(MIC_INTRINSICS)

   void SlurpIn(const char *arr, __m512i& vi, const int N_proc = N)
   {
      //_mm512_prefetch_i32gather_ps(vi, arr, 1, _MM_HINT_T0);

      const __m512    src = { 0 };
      const __mmask16 k = N_proc == N ? -1 : (1 << N_proc) - 1;

      for (int i = 0; i < kSize; ++i, arr += sizeof(T))
      {
         //_mm512_prefetch_i32gather_ps(vi, arr+2, 1, _MM_HINT_NTA);

         __m512 reg = _mm512_mask_i32gather_ps(src, k, vi, arr, 1);
         _mm512_mask_store_ps(&fArray[i*N], k, reg);
      }
   }

   /*
   // Experimental methods, SlurpIn() seems to be at least as fast.
   // See comments in mkFit/MkFitter.cc MkFitter::AddBestHit().
   void ChewIn(const char *arr, int off, int vi[N], const char *tmp,  __m512i& ui)
   {
      // This is a hack ... we know sizeof(Hit) = 64 = cache line = vector width.

      for (int i = 0; i < N; ++i)
      {
         __m512 reg = _mm512_load_ps(arr + vi[i]);
         _mm512_store_ps((void*) (tmp + 64*i), reg);
      }

      for (int i = 0; i < kSize; ++i)
      {
         __m512 reg = _mm512_i32gather_ps(ui, tmp + off + i*sizeof(T), 1);
         _mm512_store_ps(&fArray[i*N], reg);
      }
   }

   void Contaginate(const char *arr, int vi[N], const char *tmp)
   {
      // This is a hack ... we know sizeof(Hit) = 64 = cache line = vector width.

      for (int i = 0; i < N; ++i)
      {
         __m512 reg = _mm512_load_ps(arr + vi[i]);
         _mm512_store_ps((void*) (tmp + 64*i), reg);
      }
   }

   void Plexify(const char *tmp, __m512i& ui)
   {
      for (int i = 0; i < kSize; ++i)
      {
         __m512 reg = _mm512_i32gather_ps(ui, tmp + i*sizeof(T), 1);
         _mm512_store_ps(&fArray[i*N], reg);
      }
   }
   */

#else

   void SlurpIn(const char *arr, int vi[N], const int N_proc = N)
   {
      // Separate N_proc == N case (gains about 7% in fit test).
      if (N_proc == N)
      {
         for (int i = 0; i < kSize; ++i)
         {
            // Next loop vectorizes with "#pragma ivdep", but it runs slower
            // #pragma ivdep
            for (int j = 0; j < N; ++j)
            {
               fArray[i*N + j] = * (const T*) (arr + i*sizeof(T) + vi[j]);
            }
         }
      }
      else
      {
         for (int i = 0; i < kSize; ++i)
         {
            for (int j = 0; j < N_proc; ++j)
            {
               fArray[i*N + j] = * (const T*) (arr + i*sizeof(T) + vi[j]);
            }
         }
      }
   }

#endif
   
   void CopyOut(idx_t n, T *arr) const
   {
      for (idx_t i = n; i < kTotSize; i += N)
      {
         *(arr++) = fArray[i];
      }
   }
};


template<typename T, idx_t D1, idx_t D2, idx_t N> using MPlex = Matriplex<T, D1, D2, N>;


//==============================================================================
// Multiplications
//==============================================================================

template<typename T, idx_t D1, idx_t D2, idx_t D3, idx_t N>
void MultiplyGeneral(const MPlex<T, D1, D2, N>& A,
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

//------------------------------------------------------------------------------

template<typename T, idx_t D, idx_t N>
struct MultiplyCls
{
   static void Multiply(const MPlex<T, D, D, N>& A,
                        const MPlex<T, D, D, N>& B,
                        MPlex<T, D, D, N>& C)
   {
      throw std::runtime_error("general multiplication not supported, well, call MultiplyGeneral()");
   }
};

template<typename T, idx_t N>
struct MultiplyCls<T, 3, N>
{
   static void Multiply(const MPlex<T, 3, 3, N>& A,
                        const MPlex<T, 3, 3, N>& B,
                        MPlex<T, 3, 3, N>& C)
{
   const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
   const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
         T *c = C.fArray; ASSUME_ALIGNED(c, 64);

   #include "Matrix33x33.ah"
}
};

template<typename T, idx_t N>
struct MultiplyCls<T, 6, N>
{
   static void Multiply(const MPlex<T, 6, 6, N>& A,
                        const MPlex<T, 6, 6, N>& B,
                        MPlex<T, 6, 6, N>& C)
{
   const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
   const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
         T *c = C.fArray; ASSUME_ALIGNED(c, 64);

   #include "Matrix66x66.ah"
}
};

template<typename T, idx_t D, idx_t N>
void Multiply(const MPlex<T, D, D, N>& A,
              const MPlex<T, D, D, N>& B,
                    MPlex<T, D, D, N>& C)
{
   // printf("Multipl %d %d\n", D, N);

   MultiplyCls<T, D, N>::Multiply(A, B, C);
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

      T *a = A.fArray; ASSUME_ALIGNED(a, 64);

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         const TT det = a[0*N+n] * a[3*N+n] - a[2*N+n] * a[1*N+n];

         //if (determ)
         //determ[n] = det;

         const TT s   = TT(1) / det;
         const TT tmp = s * a[3*N + n];
         a[1*N+n] *= -s;
         a[2*N+n] *= -s;
         a[3*N+n]  = s * a[0*N+n];
         a[0*N+n]  = tmp;
      }
   }
};

template<typename T, idx_t N>
struct CramerInverter<T, 3, N>
{
   static void Invert(MPlex<T, 3, 3, N>& A, double *determ=0)
   {
      typedef T TT;

      T *a = A.fArray; ASSUME_ALIGNED(a, 64);

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         const TT c00 = a[4*N+n] * a[8*N+n] - a[5*N+n] * a[7*N+n];
         const TT c01 = a[5*N+n] * a[6*N+n] - a[3*N+n] * a[8*N+n];
         const TT c02 = a[3*N+n] * a[7*N+n] - a[4*N+n] * a[6*N+n];
         const TT c10 = a[7*N+n] * a[2*N+n] - a[8*N+n] * a[1*N+n];
         const TT c11 = a[8*N+n] * a[0*N+n] - a[6*N+n] * a[2*N+n];
         const TT c12 = a[6*N+n] * a[1*N+n] - a[7*N+n] * a[0*N+n];
         const TT c20 = a[1*N+n] * a[5*N+n] - a[2*N+n] * a[4*N+n];
         const TT c21 = a[2*N+n] * a[3*N+n] - a[0*N+n] * a[5*N+n];
         const TT c22 = a[0*N+n] * a[4*N+n] - a[1*N+n] * a[3*N+n];

         const TT det = a[0*N+n] * c00 + a[1*N+n] * c01 + a[2*N+n] * c02;

         //if (determ)
         //  *determ[n] = det;

         const TT s = TT(1) / det;

         a[0*N+n] = s*c00;
         a[1*N+n] = s*c10;
         a[2*N+n] = s*c20;
         a[3*N+n] = s*c01;
         a[4*N+n] = s*c11;
         a[5*N+n] = s*c21;
         a[6*N+n] = s*c02;
         a[7*N+n] = s*c12;
         a[8*N+n] = s*c22;
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

      T *a = A.fArray; ASSUME_ALIGNED(a, 64);

#pragma simd
      for (idx_t n = 0; n < N; ++n)
      {
         TT l0 = std::sqrt(T(1) / a[0*N+n]);
         TT l1 = a[3*N+n] * l0;
         TT l2 = a[4*N+n] - l1 * l1;
         l2 = std::sqrt(T(1) / l2);
         TT l3 = a[6*N+n] * l0;
         TT l4 = (a[7*N+n] - l1 * l3) * l2;
         TT l5 = a[8*N+n] - (l3 * l3 + l4 * l4);
         l5 = std::sqrt(T(1) / l5);

         // decomposition done

         l3 = (l1 * l4 * l2 - l3) * l0 * l5;
         l1 = -l1 * l0 * l2;
         l4 = -l4 * l2 * l5;

         a[0*N+n] = l3*l3 + l1*l1 + l0*l0;
         a[1*N+n] = a[3*N+n] = l3*l4 + l1*l2;
         a[4*N+n] = l4*l4 + l2*l2;
         a[2*N+n] = a[6*N+n] = l3*l5;
         a[5*N+n] = a[7*N+n] = l4*l5;
         a[8*N+n] = l5*l5;

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
