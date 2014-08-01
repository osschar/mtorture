#include "MkFitter.h"

#include "Propagation.h"

void MkFitter::SetTracksAndHits(std::vector<Track>& tracks, int beg, int end)
{
  // Assign - RAW copy -- needs to be chunked somehow

  // This is not true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Track &trk = tracks[i];

    Err[iC].Assign(itrack, trk.errors().Array());
    Par[iC].Assign(itrack, trk.parameters().Array());

    Chg(itrack, 0, 0) = trk.charge();

    for (int hi = 0; hi < Nhits; ++hi)
    {
      Hit &hit = trk.hitsVector()[hi];

      msErr[hi].Assign(itrack, hit.error().Array());
      msPar[hi].Assign(itrack, hit.parameters().Array());
    }
  }
}

void MkFitter::Fit()
{

  // Fitting loop - RAW copy -- needs to be chunked

  for (int hi = 0; hi < Nhits; ++hi)
  {
    // XXXX Note, charge is not passed (line propagation). Could be part of ctxt, too.

    propagateLineToRMPlex(Err[iC], Par[iC], msErr[hi], msPar[hi],
                          Err[iP], Par[iP],
                          updateCtx);

    updateParametersMPlex(Err[iP], Par[iP], msErr[hi], msPar[hi],
                          Err[iC], Par[iC],
                          updateCtx);
  }
}
