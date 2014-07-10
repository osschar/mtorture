#ifndef ArrayTest_h
#define ArrayTest_h

class ArrayTest
{
  float **fA;
  int     fN;

public:
  ArrayTest(int n_array, int size);

  ~ArrayTest();

  int get_n() const { return fN; }

  int copy(int n);

  int sum2(int n);
  int sum2_sqr(int n);
  int sum2_cube(int n);
  int sum2_quad(int n);
  int sum2_quint(int n);

  int sum3(int n);
  int sum3_sqr(int n);
  int sum3_cube(int n);

  // Store result back in A
  int sum2_cube_sa(int n);
  int sum2_quint_sa(int n);
  int sum3_cube_sa(int n);

};

#endif
