#include "immintrin.h"

#include <math.h>
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
#define VG(x, a, s)   _mm512_i32gather_ps(x, a, s)
#define VS(b, x, a, s)      _mm512_i32scatter_ps(b, x, a, s)
#define VSM(b, k, x, a, s)  _mm512_mask_i32scatter_ps(b, k, x, a, s)

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

  printf("Masks are not needed, but here are some facts anyway:\n");
  printf("This is how to do a vscatter of an array of length %d...\n",MM);
  const int MU = MM/NN;
  const int MO = MM%NN;
  printf("Do %d unmasked vscatter(s) of 16, then 1 masked by %d ones\n",
         MU, MO);
  int value = sizeof(unsigned short int);
  printf("Mask requires 2 bytes; unsigned short int has %d bytes\n",value);
  __mmask16 kp = (unsigned short int) (pow(2,MO) - 1);
  printf("Mask computed through int math is 0x%04x\n",kp);
  __mmask16 k = 0xffff >> (NN - MO);
  printf("Mask computed using bit shifts is 0x%04x\n",k);


  int n[NN] __attribute__((aligned(64)));
  float *pt = p;
  int j = 0;
  int m = 0;
  int s = 4;

  for (int iv = 0; iv < MM*NN; iv += NN)
  {
    __m512 r = LD(pt, 0);

    for (int i = 0; i < NN; ++i)
    {
      n[i] = j;
      j += NN;
      if (j >= NN*MM) j = ++m;
    }

    __m512i x = LDI(n, 0);
    VS(q, x, r, s);

    pt += NN;
  }

  printf(" i    p         q\n");
  for (int i = 0; i < NN; ++i)
  {
    printf("%2d %4.0f %4.0f %4.0f %4.0f %4.0f\n",
            i, p[i], p[i+NN], q[i], q[i+NN], q[i+2*NN]);
  }

  _mm_free(p);
  _mm_free(q);

  return 0;
}
