#ifndef KalmanOps_H
#define KalmanOps_H

#include "Matrix.h"

//------------------------------------------------------------------------------

namespace Matriplex
{

inline
void MultResidualsAdd(const MPlexLH& A,
                      const MPlexLV& B,
                      const MPlexHV& C,
                            MPlexLV& D)
{
   // outPar = psPar + kalmanGain*(msPar-psPar)
   //   D    =   B         A         C  -  B
   // where right half of kalman gain is 0 

   // XXX Regenerate with a script.

   typedef float T;
   const idx_t N = NN;

   const T *a = A.fArray; __assume_aligned(a, 64);
   const T *b = B.fArray; __assume_aligned(b, 64);
   const T *c = C.fArray; __assume_aligned(c, 64);
         T *d = D.fArray; __assume_aligned(d, 64);

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      // manually subrtact into local vars -- 3 of them
      float x0 = c[0 * N + n] - b[0 * N + n];
      float x1 = c[1 * N + n] - b[1 * N + n];
      float x2 = c[2 * N + n] - b[2 * N + n];

      // generate loop (can also write it manually this time, it's not much)
      d[0 * N + n] = b[0 * N + n] + a[ 0 * N + n] * x0 + a[ 1 * N + n] * x1 + a[ 2 * N + n] * x2;
      d[1 * N + n] = b[1 * N + n] + a[ 3 * N + n] * x0 + a[ 4 * N + n] * x1 + a[ 5 * N + n] * x2;
      d[2 * N + n] = b[2 * N + n] + a[ 6 * N + n] * x0 + a[ 7 * N + n] * x1 + a[ 8 * N + n] * x2;
      d[3 * N + n] = b[3 * N + n] + a[ 9 * N + n] * x0 + a[10 * N + n] * x1 + a[11 * N + n] * x2;
      d[4 * N + n] = b[4 * N + n] + a[12 * N + n] * x0 + a[13 * N + n] * x1 + a[14 * N + n] * x2;
      d[5 * N + n] = b[5 * N + n] + a[15 * N + n] * x0 + a[16 * N + n] * x1 + a[17 * N + n] * x2;
   }
}

//------------------------------------------------------------------------------

inline
void AddIntoUpperLeft3x3(const MPlexLS& A, const MPlexHS& B, MPlexHS& C)
{
   // The rest of matrix is left untouched.

   typedef float T;
   const idx_t N = NN;

   const T *a = A.fArray; __assume_aligned(a, 64);
   const T *b = B.fArray; __assume_aligned(b, 64);
         T *c = C.fArray; __assume_aligned(c, 64);

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      c[0*N+n] = a[0*N+n] + b[0*N+n];
      c[1*N+n] = a[1*N+n] + b[1*N+n];
      c[2*N+n] = a[2*N+n] + b[2*N+n];
      c[3*N+n] = a[3*N+n] + b[3*N+n];
      c[4*N+n] = a[4*N+n] + b[4*N+n];
      c[5*N+n] = a[5*N+n] + b[5*N+n];
   }
}


//------------------------------------------------------------------------------

}

#endif
