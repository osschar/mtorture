#ifndef common_h
#define common_h

#include <functional>
#include <stdexcept>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <mm_malloc.h>

typedef long long   long64;
typedef int         idx_t;

typedef std::function<long64 (int)> Func_t;

constexpr int ALIGN = 64;

template <typename X>
X* new_aligned(int n, int align=ALIGN)
{
  return (X*) _mm_malloc(sizeof(X) * n, ALIGN);
}

template <typename X>
void free_aligned(X* x)
{
  _mm_free(x);
}

#ifdef __INTEL_COMPILER
  #define ASSUME_ALIGNED(a, b) __assume_aligned(a, b)
  #define ASSUME(a, b) __assume(a == b)
#else
  #define ASSUME_ALIGNED(a, b) a = static_cast<decltype(a)>(__builtin_assume_aligned(a, b))
  #define ASSUME(a, b) __builtin_expect(a, b)
#endif

  
//------------------------------------------------------------------------------
// Globals and environment overrides
//------------------------------------------------------------------------------

template <typename X>
X get_env(const char* name, X def)
{
  char *env = getenv(name);
  return (env == 0) ? def : atof(env);
}

const double g_test_duration  = get_env("TEST_DURATION", 1.0);
const double g_pre_test_frac  = get_env("PRE_TEST_FRAC", 0.01);

const int    g_n_vec_min  = get_env("N_VEC_MIN", 1);
const int    g_n_vec_max  = get_env("N_VEC_MAX", 64 * 1024 * 1024);

#endif
