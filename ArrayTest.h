#ifndef ArrayTest_h
#define ArrayTest_h

class ArrayTest
{
  float **fA;
  int     fN;

public:
  ArrayTest(int n_array, int size);

  ~ArrayTest();

  int copy(int n);

  int sum2(int n);
  int sum2_sqr(int n);
  int sum2_cube(int n);
  int sum2_quad(int n);
  int sum2_quint(int n);

  int sum3(int n);
  int sum3_sqr(int n);
  int sum3_cube(int n);
};

#endif
