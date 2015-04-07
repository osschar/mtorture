#include "MkFitter.h"

#include "Propagation.h"

namespace
{
   inline float hipo(float x, float y) { return sqrt(x*x + y*y); }
}

void MkFitter::CheckAlignment()
{
  printf("MkFitter alignment check:\n");
  Matriplex::align_check("  Err[0]   =", &Err[0].fArray[0]);
  Matriplex::align_check("  Err[1]   =", &Err[1].fArray[0]);
  Matriplex::align_check("  Par[0]   =", &Par[0].fArray[0]);
  Matriplex::align_check("  Par[1]   =", &Par[1].fArray[0]);
  Matriplex::align_check("  msErr[0] =", &msErr[0].fArray[0]);
  Matriplex::align_check("  msPar[0] =", &msPar[0].fArray[0]);
}

void MkFitter::PrintPt(int idx)
{
  for (int i = 0; i < NN; ++i)
  {
    printf("%5.2f  ", hipo(Par[idx].At(i, 3, 0), Par[idx].At(i, 4, 0)));
  }
}

//==============================================================================

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

//==============================================================================

void MkFitter::InputIntrTracksAndHits(std::vector<Track>& tracks, int beg, int end)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Track &trk = tracks[i];

    Err[iC].CopyInIntr(itrack, trk.errors().Array());
    Par[iC].CopyInIntr(itrack, trk.parameters().Array());

    Chg(itrack, 0, 0) = trk.charge();

    for (int hi = 0; hi < Nhits; ++hi)
    {
      Hit &hit = trk.hitsVector()[hi];

      msErr[hi].CopyInIntr(itrack, hit.error().Array());
      msPar[hi].CopyInIntr(itrack, hit.parameters().Array());
    }
  }
}

//==============================================================================

void MkFitter::InputContigTracksAndHits(std::vector<Track>& tracks, int beg, int end)
{
  // Assign track parameters to initial state and copy hit values in -
  // BUT - store items belonging to the same track in contiguous memory,
  // i.e., do not put the Matriplex constituents into Matriplex order.
  // This method differs from InputTracksAndHits in that it requires a
  // subsequent call to PlexifyTracksAndHits to create a valid MkFitter.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Track &trk = tracks[i];

    Err[iC].CopyInContig(itrack, trk.errors().Array());
    Par[iC].CopyInContig(itrack, trk.parameters().Array());

    Chg(itrack, 0, 0) = trk.charge();

    for (int hi = 0; hi < Nhits; ++hi)
    {
      Hit &hit = trk.hitsVector()[hi];

      msErr[hi].CopyInContig(itrack, hit.error().Array());
      msPar[hi].CopyInContig(itrack, hit.parameters().Array());
    }
  }
}

//==============================================================================

void MkFitter::PlexifyTracksAndHits(MkFitter& MkFContig)
{
  Err[iC].Plexify(MkFContig.Err[iC].fArray);
  Par[iC].Plexify(MkFContig.Par[iC].fArray);

  for (int i = 0; i < NN; ++i)
  {
    Chg(i, 0, 0) = MkFContig.Chg(i, 0, 0);
  }

  for (int hi = 0; hi < Nhits; ++hi)
  {
    msErr[hi].Plexify(MkFContig.msErr[hi].fArray);
    msPar[hi].Plexify(MkFContig.msPar[hi].fArray);
  }
}

//==============================================================================

void MkFitter::PlexifyIntrTracksAndHits(MkFitter& MkFContig)
{
  Err[iC].PlexifyIntr(MkFContig.Err[iC].fArray);
  Par[iC].PlexifyIntr(MkFContig.Par[iC].fArray);

  for (int i = 0; i < NN; ++i)
  {
    Chg(i, 0, 0) = MkFContig.Chg(i, 0, 0);
  }

  for (int hi = 0; hi < Nhits; ++hi)
  {
    msErr[hi].PlexifyIntr(MkFContig.msErr[hi].fArray);
    msPar[hi].PlexifyIntr(MkFContig.msPar[hi].fArray);
  }
}

//==============================================================================

void MkFitter::PlexifyIntr2TracksAndHits(MkFitter& MkFContig)
{
  Err[iC].PlexifyIntr2(MkFContig.Err[iC].fArray);
  Par[iC].PlexifyIntr2(MkFContig.Par[iC].fArray);

  for (int i = 0; i < NN; ++i)
  {
    Chg(i, 0, 0) = MkFContig.Chg(i, 0, 0);
  }

  for (int hi = 0; hi < Nhits; ++hi)
  {
    msErr[hi].PlexifyIntr2(MkFContig.msErr[hi].fArray);
    msPar[hi].PlexifyIntr2(MkFContig.msPar[hi].fArray);
  }
}

//==============================================================================

#ifdef MKLOPT
void MkFitter::PlexifyMKLOutTracksAndHits(MkFitter& MkFContig)
{
  Err[iC].PlexifyMKLOut(MkFContig.Err[iC].fArray);
  Par[iC].PlexifyMKLOut(MkFContig.Par[iC].fArray);

  for (int i = 0; i < NN; ++i)
  {
    Chg(i, 0, 0) = MkFContig.Chg(i, 0, 0);
  }

  for (int hi = 0; hi < Nhits; ++hi)
  {
    msErr[hi].PlexifyMKLOut(MkFContig.msErr[hi].fArray);
    msPar[hi].PlexifyMKLOut(MkFContig.msPar[hi].fArray);
  }
}
#endif

//==============================================================================

void MkFitter::InputTracksOnly(std::vector<Track>& tracks, int beg, int end)
{
  // Assign track parameters to initial state, do NOT copy hit values.
  // Used for benchmarking the fitting with less "copy-in" load.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Track &trk = tracks[i];

    Err[iC].CopyIn(itrack, trk.errors().Array());
    Par[iC].CopyIn(itrack, trk.parameters().Array());

    Chg(itrack, 0, 0) = trk.charge();
  }
}

void MkFitter::FitTracks()
{
  // Fitting loop.

  for (int hi = 0; hi < Nhits; ++hi)
  {
    // Note, charge is not passed (line propagation).
    // propagateLineToRMPlex(Err[iC], Par[iC], msErr[hi], msPar[hi],
    //                       Err[iP], Par[iP]);

    propagateHelixToRMPlex(Err[iC], Par[iC], Chg, msPar[hi],
                           Err[iP], Par[iP]);

    updateParametersMPlex(Err[iP], Par[iP], msErr[hi], msPar[hi],
                          Err[iC], Par[iC]);
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
    Err[iC].CopyOut(itrack, tracks[i].errors().Array());
    Par[iC].CopyOut(itrack, tracks[i].parameters().Array());

    tracks[i].charge() = Chg(itrack, 0, 0);

    // XXXXX chi2 is not set (also not in SMatrix fit, it seems)
  }
}
