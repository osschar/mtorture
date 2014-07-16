# include "MPlexTest.h"

const int MPlexTest::sNMul = MPT_DIM*MPT_DIM*(2*MPT_DIM - 1) * MPT_SIZE;

//==============================================================================

MPlexTest::MPlexTest(int n_array, int size)
{
  fMPV = new MPV*[n_array];
  fN   = n_array;

  for (int i = 0; i < fN; ++i)
  {
    fMPV[i] = new MPV(size);
  }
}

MPlexTest::~MPlexTest()
{
  for (int i = 0; i < fN; ++i)
  {
    delete fMPV[i];
  }
  delete fMPV;
}

//==============================================================================

int MPlexTest::mult2(int n_vec)
{
  Matriplex::MultiplyUnrolled(*fMPV[0], *fMPV[1], *fMPV[2], n_vec);

  return sNMul * n_vec;
}

int MPlexTest::mult2_3out(int n_vec)
{
  Matriplex::Multiply(*fMPV[0], *fMPV[1], *fMPV[2], n_vec);
  Matriplex::Multiply(*fMPV[1], *fMPV[2], *fMPV[0], n_vec);
  Matriplex::Multiply(*fMPV[2], *fMPV[0], *fMPV[1], n_vec);

  return 3 * sNMul * n_vec;
}

int MPlexTest::mult2_3in(int n_vec)
{
  Matriplex::Multiply3in(*fMPV[0], *fMPV[1], *fMPV[2], n_vec);

  return 3 * sNMul * n_vec;
}

//==============================================================================

int MPlexTest::inv_cramer(int n_vec)
{
  assert(MPT_DIM == 3);

  Matriplex::InvertCramer(*fMPV[0], n_vec);

  // 3x3: 27 + 5 + 1 + 9 = 42 
  return 42 * MPT_SIZE * n_vec;
}

int MPlexTest::inv_cholesky(int n_vec)
{
  assert(MPT_DIM == 3);

  Matriplex::InvertCholesky(*fMPV[0], n_vec);

  // 3x3: 32 ops + 3 sqrts
  // sqrt: latency 8, throughput 2 (have to wait on each of them!)
  // what do i do then? count each as 4 (every other clock!).
  return 44 * MPT_SIZE * n_vec;
}
