#include "MatriplexCommon.h"

#include "fittest.h"

#include "MkFitter.h"

#include "Timing.h"

#include <limits>

std::vector<Track> simtracks;

std::vector<Track> smat_tracks;
std::vector<Track> plex_tracks;

//==============================================================================

MkFitter *mkfp;

const int Nhits = 10; // XXXXX ARGH !!!! What if there's a missing / double layer?

const int Nloop = 1;

long64 single_run(int n_tracks)
{
  int theEnd = n_tracks;

  for (int itrack = 0; itrack < theEnd; itrack += NN)
  {
    int end = std::min(itrack + NN, theEnd);

    mkfp->InputTracksAndHits(simtracks, itrack, end);

    for (int x = 0; x < Nloop; ++x)
    {
      if (x != 0)
      {
        mkfp->InputTracksOnly(simtracks, itrack, end);
      }

      mkfp->FitTracks();
    }

#ifndef NO_ROOT
    mkfp->OutputFittedTracks(rectracks, itrack, end);
#endif
  }

  return long64(Nloop) * n_tracks * 374 * Nhits;
}

void test_matriplex()
{
  int Nmin  = 16;
  int Nmax  = 32*1024; //32 * 1024;

  generateTracks(simtracks, Nmax);

  mkfp = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Nhits);

  Timing t([&](int n_vec)
           {
             return single_run(n_vec);
           });

  t.print_tuple_header();

  // Warm up
  t.time_loop(4096, 4);

  for (int n = Nmin; n <= Nmax; n *= 2)
  {
    // Nloop = std::numeric_limits<int>::max() / (374 * Nhits * n);

    // printf("XXX n=%d, nloop=%d\n", n, Nloop);

    // t.auto_time_loop(n, 2);
    
    t.time_loop(n, 1);

    t.print(n);
  }

  _mm_free(mkfp);
}

//------------------------------------------------------------------------------

void test_standard()
{
  int  Ntracks  = 16 * 1024;
  bool saveTree = false;

  generateTracks(simtracks, Ntracks);

  smat_tracks.reserve(simtracks.size());
  plex_tracks.resize (simtracks.size());


  double tmp, tsm;

  tsm = runFittingTest(simtracks, smat_tracks);

  tmp = runFittingTestPlex(simtracks, plex_tracks);

  printf("SMatrix = %.3f   Matriplex = %.3f   ---   SM/MP = %.3f\n", tsm, tmp, tsm / tmp);

#ifndef NO_ROOT
  make_validation_tree("validation-smat.root", simtracks, smat_tracks);
  make_validation_tree("validation-plex.root", simtracks, plex_tracks);
#endif
}

//==============================================================================

int main()
{
  test_standard();

  // test_matriplex();

  return 0;
}
