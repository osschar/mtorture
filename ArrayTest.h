#ifndef ArrayTest_h
#define ArrayTest_h

#include "common.h"

class ArrayTest
{
  float **fA;
  int     fN;

public:
  ArrayTest(int n_array, int size);

  ~ArrayTest();

  int get_n() const { return fN; }

  long64 copy(int n);

  long64 sum2(int n);
  long64 sum2_sqr(int n);
  long64 sum2_cube(int n);
  long64 sum2_quad(int n);
  long64 sum2_quint(int n);

  long64 sum3(int n);
  long64 sum3_sqr(int n);
  long64 sum3_cube(int n);

  // Store result back in A
  long64 sum2_cube_sa(int n);
  long64 sum2_quint_sa(int n);
  long64 sum3_cube_sa(int n);

};

#endif
