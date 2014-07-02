#ifndef common_h
#define common_h

#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>

typedef long long   long64;
typedef int         idx_t;

const int ALIGN   = 64;

// Set this to 8 for AVX, 16 for MIC
const idx_t Sfac = 1;
#ifdef __MIC__

const idx_t S = 16 * Sfac;
#else
const idx_t S = 8  * Sfac;
#endif


template <typename X>
X* new_sth(int n, int align=ALIGN)
{
  return (X*) _mm_malloc(sizeof(X) * n, ALIGN);
}

template <typename X>
void free_sth(X* x)
{
  _mm_free(x);
}

#endif
