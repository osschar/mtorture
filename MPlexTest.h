#ifndef MPlexText_h
#define MPlexText_h

#include "Matriplex/Matriplex.h"
#include "Matriplex/MatriplexVector.h"

#ifndef MPT_DIM
#define MPT_DIM 6
#endif
#ifndef MPT_SIZE
#define MPT_SIZE 16
#endif

class MPlexTest
{
  typedef Matriplex::Matriplex<float, MPT_DIM, MPT_DIM, MPT_SIZE> MP;
  typedef Matriplex::MatriplexVector<MP>                          MPV;

  MPV **fMPV;
  int   fN;

  static const int sNMul;

public:
  MPlexTest(int n_array, int size);
  ~MPlexTest();

  int mult2(int n_vec);
  int mult2_3out(int n_vec);
  int mult2_3in (int n_vec);
};

#endif;
