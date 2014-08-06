#include "MkFitter.h"

#include "Propagation.h"

void MkFitter::InputTracksAndHits(std::vector<Track>& tracks, int beg, int end)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Track &trk = tracks[i];

    Err[iC].CopyIn(itrack, trk.errors().Array());
    Par[iC].CopyIn(itrack, trk.parameters().Array());

    Chg(itrack, 0, 0) = trk.charge();

    for (int hi = 0; hi < Nhits; ++hi)
    {
      Hit &hit = trk.hitsVector()[hi];

      msErr[hi].CopyIn(itrack, hit.error().Array());
      msPar[hi].CopyIn(itrack, hit.parameters().Array());
    }
  }
}

void MkFitter::FitTracks()
{
  // Fitting loop.

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

  // XXXXX What's with chi2?
}

void MkFitter::OutputFittedTracks(std::vector<Track>& tracks, int beg, int end)
{
  // Copies last track parameters (updated) into Track objects.
  // The tracks vector should be resized to allow direct copying.

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Err[iP].CopyOut(itrack, tracks[i].errors().Array());
    Par[iP].CopyOut(itrack, tracks[i].parameters().Array());

    tracks[i].charge() = Chg(itrack, 0, 0);

    // XXXXX chi2 is not set (also not in SMatrix fit, it seems)
  }
}
