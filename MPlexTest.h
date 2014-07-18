#ifndef MPlexText_h
#define MPlexText_h

#include "Matriplex/Matriplex.h"
#include "Matriplex/MatriplexSym.h"
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

  typedef Matriplex::MatriplexSym<float, MPT_DIM, MPT_SIZE>       MPS;
  typedef Matriplex::MatriplexVector<MPS>                         MPVS;

  MPV  **fMPV;
  MPVS **fMPVS;
  int   fN;
  int   fNS;

  static const int sNMul;
  // static const int sNMulSym;

public:
  MPlexTest(int n_mp, int n_mp_sym, int size);
  ~MPlexTest();

  int mult2(int n_vec);
  int mult2_3out(int n_vec);
  int mult2_3in (int n_vec);

  int inv_cramer(int n_vec);
  int inv_cholesky(int n_vec);

  int mult2_sym(int n_vec);

  int inv_cramer_sym(int n_vec);
  int inv_cholesky_sym(int n_vec);
};

#endif;
