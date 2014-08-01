#ifndef MkFitter_h
#define MkFitter_h

#include "Matrix.h"
#include "KalmanUtils.h"

class MkFitter
{
  MPlexSS psErr;  MPlexMV psPar;
  MPlexSS outErr; MPlexMV outPar;

  MPlexSS Err[2];
  MPlexMV Par[2];

  MPlexQI Chg;

  int iC = 0; // current
  int iP = 1; // propagated


  std::vector<MPlexSS> msErr;
  std::vector<MPlexMV> msPar;

  updateParametersContext updateCtx;

  int Nhits;

public:
  MkFitter(int n_hits) : Nhits(n_hits)
  {
    msErr.resize(Nhits);
    msPar.resize(Nhits);
  }

  void SetTracksAndHits(std::vector<Track>& tracks, int beg, int end);
  void Fit();
};

#endif
