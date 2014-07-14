#ifndef MPlexText_h
#define MPlexText_h

#include "Matriplex/Matriplex.h"
#include "Matriplex/MatriplexVector.h"

class MPlexTest
{
  typedef Matriplex::Matriplex<float, 6, 6, 16>       MP;// XXXX
  typedef Matriplex::MatriplexVector<float, 6, 6, 16> MPV;// XXXX

  MPV **fMPV;
  int   fN;

public:
  MPlexTest(int n_array, int size);
  ~MPlexTest();

  int mult2(int n_vec);
  int mult2_3out(int n_vec);
  int mult2_3in (int n_vec);
};

#endif;
