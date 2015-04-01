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
#define VS(b, x, a, s)   _mm512_i32scatter_ps(b, x, a, s);


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

  int s = 1;

  float *pt = p;

  for (int n = 0; n < MM*NN; n += MM)
  {
    __m512 r = LD(pt, 0);

    int c[NN];
    // define position n within each vector in Matriplex
    for (int i = 0; i < MM; ++i)
    {
       c[i] = n + i*NN;
    }
    printf("before load %d\n",n);
    __m512i x = LDI(c, 0);
    printf("after load %d\n",n);
    printf("before scatter %d\n",n);
    VS(q, x, r, s);
    printf("after scatter %d\n",n);
    pt += MM;
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
