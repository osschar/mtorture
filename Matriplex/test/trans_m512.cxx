#include "immintrin.h"

#include <cstdio>

const int MM = 21;
const int NN = 16;

#define LD(a, i)      _mm512_load_ps(&a[i*16])
#define ADD(a, b)     _mm512_add_ps(a, b) 
#define MUL(a, b)     _mm512_mul_ps(a, b)
#define FMA(a, b, v)  _mm512_fmadd_ps(a, b, v)
#define ST(a, i, r)   _mm512_store_ps(&a[i*16], r)
#define LDI(n, i)     _mm512_load_epi32(&n[i*16])
#define ADDI(n, m)    _mm512_add_epi32(n, m) 
#define VG(x, a, s)   _mm512_i32gather_ps(x, a, s);

// Can even be global!
__m512i all_uno = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

int main()
{
  float *p = (float*) _mm_malloc(MM*NN*sizeof(float), 64);
  float *q = (float*) _mm_malloc(MM*NN*sizeof(float), 64);

  for (int i = 0; i < MM*NN; ++i)
  {
    p[i] = i;
  }

  int n[NN];

  for (int i = 0; i < NN; ++i)  // NN starts of matrices
  {
    n[i] = i * MM;
  }
  __m512i xn = LDI(n, 0);
  int s = 4;

  float *pt = p, *qt = q;
  for (int i = 0; i < MM*NN; i += NN)
  {
    __m512 c = VG(xn, pt, s);
    ST(qt, 0, c);
    ++pt;
    qt += NN;
  }

  for (int i = 0; i < 16; ++i)
  {
    printf("%2d %4.0f %4.0f %4.0f %4.0f %4.0f\n",
            i, p[i], p[i+NN], q[i], q[i+NN], q[i+2*NN]);
  }

  _mm_free(p);
  _mm_free(q);

  return 0;
}
