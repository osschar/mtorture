#ifndef MPlexText_h
#define MPlexText_h

#include "common.h"

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

  // ----------------------------------------------------------------

  MPV&  MPlexVec   (int i) { return *fMPV[i];  }
  MPVS& MPlexVecSym(int i) { return *fMPVS[i]; }

  MP&   MPlex   (int i, int m) { return (*fMPV[i])[m];  }
  MPS&  MPlexSym(int i, int m) { return (*fMPVS[i])[m]; }

  Func_t name_to_func(const std::string& name);

  // ----------------------------------------------------------------

  long64 mult2(int n_vec);
  long64 mult2_3out(int n_vec);
  long64 mult2_3in (int n_vec);

  long64 inv_cramer(int n_vec);
  long64 inv_cholesky(int n_vec);

  long64 mult2_sym(int n_vec);

  long64 inv_cramer_sym(int n_vec);
  long64 inv_cholesky_sym(int n_vec);
};

#endif
