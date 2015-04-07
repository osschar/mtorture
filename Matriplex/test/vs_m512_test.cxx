#include "immintrin.h"

#include <cstdio>

const int NN = 16;

#define LD(a, i)      _mm512_load_ps(&a[i*16])
#define ADD(a, b)     _mm512_add_ps(a, b) 
#define MUL(a, b)     _mm512_mul_ps(a, b)
#define FMA(a, b, v)  _mm512_fmadd_ps(a, b, v)
#define ST(a, i, r)   _mm512_store_ps(&a[i*16], r)
#define LDI(n, i)     _mm512_load_epi32(&n[i*16])
#define VG(x, a, s)   _mm512_i32gather_ps(x, a, s);
#define VS(b, x, a, s)      _mm512_i32scatter_ps(b, x, a, s)

// Can even be global!
__m512 all_ones = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

int main()
{
  float *p = (float*) _mm_malloc(2*NN*sizeof(float), 64);
  float *q = (float*) _mm_malloc(2*NN*sizeof(float), 64);

  int n[NN] __attribute__((aligned(64)));

  for (int i = 0; i < NN; ++i)
  {
    p[i] = i;
    n[i] = i;
  }

  //int s = 4;
  int s = 8;
  __m512i x = LDI(n, 0);

  __m512 a = LD(p, 0);
  //__m512 a = VG(x, p, s);
  __m512 b = { 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 };//LD(p, 1);

  //b = all_ones;

  __m512 c = ADD(a, b);

  //ST(q, 0, c);
  VS(q, x, c, s);

  printf(" i    p         q\n");
  for (int i = 0; i < 16; ++i)
  {
    printf("%2d %4.0f %4.0f %4.0f %4.0f\n", i, p[i], p[i+16], q[i], q[i+16]);
  }

  _mm_free(p);
  _mm_free(q);

  return 0;
}
