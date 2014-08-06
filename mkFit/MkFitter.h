#ifndef MkFitter_h
#define MkFitter_h

#include "Matrix.h"
#include "KalmanUtils.h"

#define MAX_HITS 10


class MkFitter
{
  MPlexLS Err[2];
  MPlexLV Par[2];

  MPlexQI Chg;

  updateParametersContext updateCtx;

  MPlexHS msErr[MAX_HITS];
  MPlexHV msPar[MAX_HITS];

  // Indices into Err and Par arrays.
  // Thought I'll have to flip between them ...
  const int iC = 0; // current
  const int iP = 1; // propagated

  int Nhits;

public:
  MkFitter(int n_hits) : Nhits(n_hits)
  {
    // XXXX Eventually dynamically allocate measurement arrays
    // msErr.resize(Nhits);
    // msPar.resize(Nhits);

    printf("MkFitter alignment check:\n");
    Matriplex::align_check("  Err[0]   =", &Err[0].fArray[0]);
    Matriplex::align_check("  Err[1]   =", &Err[1].fArray[0]);
    Matriplex::align_check("  Par[0]   =", &Par[0].fArray[0]);
    Matriplex::align_check("  Par[1]   =", &Par[1].fArray[0]);
    Matriplex::align_check("  msErr[0] =", &msErr[0].fArray[0]);
    Matriplex::align_check("  msPar[0] =", &msPar[0].fArray[0]);
  }

  void InputTracksAndHits(std::vector<Track>& tracks, int beg, int end);
  void FitTracks();
  void OutputFittedTracks(std::vector<Track>& tracks, int beg, int end);

} __attribute__((aligned(64)));

#endif
