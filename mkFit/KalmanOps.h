#ifndef KalmanOps_H
#define KalmanOps_H

#include "MatriplexSym.h"

//------------------------------------------------------------------------------

namespace Matriplex
{

typedef float T; // XXXXX This is horrible!!!

/*
inline
void MultForKalmanGain(const MPlexLS& A,
                       const MPlexLS& B,
                             MPlexLH& C)
{
   // calculate Kalman gain -- multiplication where B / resErr is only populated
   // in upper-left 3x3
   // kalmanGain = propErr * resErrInverse
   //     C      =    A    *      B

   // XXX Regenerate with a script.

   const idx_t N = NN;

   const T *a = A.fArray; __assume_aligned(a, 64);
   const T *b = B.fArray; __assume_aligned(b, 64);
         T *c = C.fArray; __assume_aligned(c, 64);

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      C.fArray[0 * N + n] = A.fArray[0 * N + n] * B.fArray[0 * N + n] + A.fArray[1 * N + n] * B.fArray[1 * N + n] + A.fArray[3 * N + n] * B.fArray[3 * N + n];
      C.fArray[1 * N + n] = A.fArray[0 * N + n] * B.fArray[1 * N + n] + A.fArray[1 * N + n] * B.fArray[2 * N + n] + A.fArray[3 * N + n] * B.fArray[4 * N + n];
      C.fArray[2 * N + n] = A.fArray[0 * N + n] * B.fArray[3 * N + n] + A.fArray[1 * N + n] * B.fArray[4 * N + n] + A.fArray[3 * N + n] * B.fArray[5 * N + n];
      C.fArray[6 * N + n] = A.fArray[1 * N + n] * B.fArray[0 * N + n] + A.fArray[2 * N + n] * B.fArray[1 * N + n] + A.fArray[4 * N + n] * B.fArray[3 * N + n];
      C.fArray[7 * N + n] = A.fArray[1 * N + n] * B.fArray[1 * N + n] + A.fArray[2 * N + n] * B.fArray[2 * N + n] + A.fArray[4 * N + n] * B.fArray[4 * N + n];
      C.fArray[8 * N + n] = A.fArray[1 * N + n] * B.fArray[3 * N + n] + A.fArray[2 * N + n] * B.fArray[4 * N + n] + A.fArray[4 * N + n] * B.fArray[5 * N + n];
      C.fArray[12 * N + n] = A.fArray[3 * N + n] * B.fArray[0 * N + n] + A.fArray[4 * N + n] * B.fArray[1 * N + n] + A.fArray[5 * N + n] * B.fArray[3 * N + n];
      C.fArray[13 * N + n] = A.fArray[3 * N + n] * B.fArray[1 * N + n] + A.fArray[4 * N + n] * B.fArray[2 * N + n] + A.fArray[5 * N + n] * B.fArray[4 * N + n];
      C.fArray[14 * N + n] = A.fArray[3 * N + n] * B.fArray[3 * N + n] + A.fArray[4 * N + n] * B.fArray[4 * N + n] + A.fArray[5 * N + n] * B.fArray[5 * N + n];
      C.fArray[18 * N + n] = A.fArray[6 * N + n] * B.fArray[0 * N + n] + A.fArray[7 * N + n] * B.fArray[1 * N + n] + A.fArray[8 * N + n] * B.fArray[3 * N + n];
      C.fArray[19 * N + n] = A.fArray[6 * N + n] * B.fArray[1 * N + n] + A.fArray[7 * N + n] * B.fArray[2 * N + n] + A.fArray[8 * N + n] * B.fArray[4 * N + n];
      C.fArray[20 * N + n] = A.fArray[6 * N + n] * B.fArray[3 * N + n] + A.fArray[7 * N + n] * B.fArray[4 * N + n] + A.fArray[8 * N + n] * B.fArray[5 * N + n];
      C.fArray[24 * N + n] = A.fArray[10 * N + n] * B.fArray[0 * N + n] + A.fArray[11 * N + n] * B.fArray[1 * N + n] + A.fArray[12 * N + n] * B.fArray[3 * N + n];
      C.fArray[25 * N + n] = A.fArray[10 * N + n] * B.fArray[1 * N + n] + A.fArray[11 * N + n] * B.fArray[2 * N + n] + A.fArray[12 * N + n] * B.fArray[4 * N + n];
      C.fArray[26 * N + n] = A.fArray[10 * N + n] * B.fArray[3 * N + n] + A.fArray[11 * N + n] * B.fArray[4 * N + n] + A.fArray[12 * N + n] * B.fArray[5 * N + n];
      C.fArray[30 * N + n] = A.fArray[15 * N + n] * B.fArray[0 * N + n] + A.fArray[16 * N + n] * B.fArray[1 * N + n] + A.fArray[17 * N + n] * B.fArray[3 * N + n];
      C.fArray[31 * N + n] = A.fArray[15 * N + n] * B.fArray[1 * N + n] + A.fArray[16 * N + n] * B.fArray[2 * N + n] + A.fArray[17 * N + n] * B.fArray[4 * N + n];
      C.fArray[32 * N + n] = A.fArray[15 * N + n] * B.fArray[3 * N + n] + A.fArray[16 * N + n] * B.fArray[4 * N + n] + A.fArray[17 * N + n] * B.fArray[5 * N + n];
   }
}
*/

//------------------------------------------------------------------------------

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

   const idx_t N = NN;

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      // manually subrtact into local vars -- 3 of them
      float d0 = C.fArray[0 * N + n] - B.fArray[0 * N + n];
      float d1 = C.fArray[1 * N + n] - B.fArray[1 * N + n];
      float d2 = C.fArray[2 * N + n] - B.fArray[2 * N + n];

      // generate loop (can also write it manually this time, it's not much)
      D.fArray[0 * N + n] = B.fArray[0 * N + n] + A.fArray[ 0 * N + n] * d0 + A.fArray[ 1 * N + n] * d1 + A.fArray[ 2 * N + n] * d2;
      D.fArray[1 * N + n] = B.fArray[1 * N + n] + A.fArray[ 3 * N + n] * d0 + A.fArray[ 4 * N + n] * d1 + A.fArray[ 5 * N + n] * d2;
      D.fArray[2 * N + n] = B.fArray[2 * N + n] + A.fArray[ 6 * N + n] * d0 + A.fArray[ 7 * N + n] * d1 + A.fArray[ 8 * N + n] * d2;
      D.fArray[3 * N + n] = B.fArray[3 * N + n] + A.fArray[ 9 * N + n] * d0 + A.fArray[10 * N + n] * d1 + A.fArray[11 * N + n] * d2;
      D.fArray[4 * N + n] = B.fArray[4 * N + n] + A.fArray[12 * N + n] * d0 + A.fArray[13 * N + n] * d1 + A.fArray[14 * N + n] * d2;
      D.fArray[5 * N + n] = B.fArray[5 * N + n] + A.fArray[15 * N + n] * d0 + A.fArray[16 * N + n] * d1 + A.fArray[17 * N + n] * d2;
   }
}

//------------------------------------------------------------------------------

inline
void FinalKalmanErr(const MPlexLS& A,
                    const MPlexLL& B,
                          MPlexLS& C)
{
   // propErr - kalmanGain*propErr
   //    A    -      B    *   A
   // where right half of B is 0

   // XXX Regenerate with a script.

   const idx_t N = NN;

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      C.fArray[0 * N + n] = A.fArray[0 * N + n] - B.fArray[0 * N + n] * A.fArray[0 * N + n] - B.fArray[1 * N + n] * A.fArray[1 * N + n] - B.fArray[2 * N + n] * A.fArray[3 * N + n];
      C.fArray[1 * N + n] = A.fArray[1 * N + n] - B.fArray[6 * N + n] * A.fArray[0 * N + n] - B.fArray[7 * N + n] * A.fArray[1 * N + n] - B.fArray[8 * N + n] * A.fArray[3 * N + n];
      C.fArray[2 * N + n] = A.fArray[2 * N + n] - B.fArray[6 * N + n] * A.fArray[1 * N + n] - B.fArray[7 * N + n] * A.fArray[2 * N + n] - B.fArray[8 * N + n] * A.fArray[4 * N + n];
      C.fArray[3 * N + n] = A.fArray[3 * N + n] - B.fArray[12 * N + n] * A.fArray[0 * N + n] - B.fArray[13 * N + n] * A.fArray[1 * N + n] - B.fArray[14 * N + n] * A.fArray[3 * N + n];
      C.fArray[4 * N + n] = A.fArray[4 * N + n] - B.fArray[12 * N + n] * A.fArray[1 * N + n] - B.fArray[13 * N + n] * A.fArray[2 * N + n] - B.fArray[14 * N + n] * A.fArray[4 * N + n];
      C.fArray[5 * N + n] = A.fArray[5 * N + n] - B.fArray[12 * N + n] * A.fArray[3 * N + n] - B.fArray[13 * N + n] * A.fArray[4 * N + n] - B.fArray[14 * N + n] * A.fArray[5 * N + n];
      C.fArray[6 * N + n] = A.fArray[6 * N + n] - B.fArray[18 * N + n] * A.fArray[0 * N + n] - B.fArray[19 * N + n] * A.fArray[1 * N + n] - B.fArray[20 * N + n] * A.fArray[3 * N + n];
      C.fArray[7 * N + n] = A.fArray[7 * N + n] - B.fArray[18 * N + n] * A.fArray[1 * N + n] - B.fArray[19 * N + n] * A.fArray[2 * N + n] - B.fArray[20 * N + n] * A.fArray[4 * N + n];
      C.fArray[8 * N + n] = A.fArray[8 * N + n] - B.fArray[18 * N + n] * A.fArray[3 * N + n] - B.fArray[19 * N + n] * A.fArray[4 * N + n] - B.fArray[20 * N + n] * A.fArray[5 * N + n];
      C.fArray[9 * N + n] = A.fArray[9 * N + n] - B.fArray[18 * N + n] * A.fArray[6 * N + n] - B.fArray[19 * N + n] * A.fArray[7 * N + n] - B.fArray[20 * N + n] * A.fArray[8 * N + n];
      C.fArray[10 * N + n] = A.fArray[10 * N + n] - B.fArray[24 * N + n] * A.fArray[0 * N + n] - B.fArray[25 * N + n] * A.fArray[1 * N + n] - B.fArray[26 * N + n] * A.fArray[3 * N + n];
      C.fArray[11 * N + n] = A.fArray[11 * N + n] - B.fArray[24 * N + n] * A.fArray[1 * N + n] - B.fArray[25 * N + n] * A.fArray[2 * N + n] - B.fArray[26 * N + n] * A.fArray[4 * N + n];
      C.fArray[12 * N + n] = A.fArray[12 * N + n] - B.fArray[24 * N + n] * A.fArray[3 * N + n] - B.fArray[25 * N + n] * A.fArray[4 * N + n] - B.fArray[26 * N + n] * A.fArray[5 * N + n];
      C.fArray[13 * N + n] = A.fArray[13 * N + n] - B.fArray[24 * N + n] * A.fArray[6 * N + n] - B.fArray[25 * N + n] * A.fArray[7 * N + n] - B.fArray[26 * N + n] * A.fArray[8 * N + n];
      C.fArray[14 * N + n] = A.fArray[14 * N + n] - B.fArray[24 * N + n] * A.fArray[10 * N + n] - B.fArray[25 * N + n] * A.fArray[11 * N + n] - B.fArray[26 * N + n] * A.fArray[12 * N + n];
      C.fArray[15 * N + n] = A.fArray[15 * N + n] - B.fArray[30 * N + n] * A.fArray[0 * N + n] - B.fArray[31 * N + n] * A.fArray[1 * N + n] - B.fArray[32 * N + n] * A.fArray[3 * N + n];
      C.fArray[16 * N + n] = A.fArray[16 * N + n] - B.fArray[30 * N + n] * A.fArray[1 * N + n] - B.fArray[31 * N + n] * A.fArray[2 * N + n] - B.fArray[32 * N + n] * A.fArray[4 * N + n];
      C.fArray[17 * N + n] = A.fArray[17 * N + n] - B.fArray[30 * N + n] * A.fArray[3 * N + n] - B.fArray[31 * N + n] * A.fArray[4 * N + n] - B.fArray[32 * N + n] * A.fArray[5 * N + n];
      C.fArray[18 * N + n] = A.fArray[18 * N + n] - B.fArray[30 * N + n] * A.fArray[6 * N + n] - B.fArray[31 * N + n] * A.fArray[7 * N + n] - B.fArray[32 * N + n] * A.fArray[8 * N + n];
      C.fArray[19 * N + n] = A.fArray[19 * N + n] - B.fArray[30 * N + n] * A.fArray[10 * N + n] - B.fArray[31 * N + n] * A.fArray[11 * N + n] - B.fArray[32 * N + n] * A.fArray[12 * N + n];
      C.fArray[20 * N + n] = A.fArray[20 * N + n] - B.fArray[30 * N + n] * A.fArray[15 * N + n] - B.fArray[31 * N + n] * A.fArray[16 * N + n] - B.fArray[32 * N + n] * A.fArray[17 * N + n];
   }
}

//------------------------------------------------------------------------------

inline
void AddIntoUpperLeft3x3(const MPlexLS& A, const MPlexHS& B, MPlexHS& C)
{
   // The rest of matrix is left untouched.

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
