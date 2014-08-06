#ifndef MkFitter_h
#define MkFitter_h

#include "Matrix.h"
#include "KalmanUtils.h"

class MkFitter
{
  MPlexLS Err[2];
  MPlexLV Par[2];

  MPlexQI Chg;

  updateParametersContext updateCtx;

  std::vector<MPlexHS> msErr;
  std::vector<MPlexHV> msPar;

  // Indices into Err and Par arrays.
  // Thought I'll have to flip between them ...
  const int iC = 0; // current
  const int iP = 1; // propagated

  int Nhits;

public:
  MkFitter(int n_hits) : Nhits(n_hits)
  {
    msErr.resize(Nhits);
    msPar.resize(Nhits);
  }

  void InputTracksAndHits(std::vector<Track>& tracks, int beg, int end);
  void FitTracks();
  void OutputFittedTracks(std::vector<Track>& tracks, int beg, int end);
};

#endif
