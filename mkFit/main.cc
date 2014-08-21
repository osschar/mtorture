#include "MatriplexCommon.h"

#include "fittest.h"

#include "MkFitter.h"

#include "Timing.h"

#include <limits>

#ifndef NUM_THREADS
#define NUM_THREADS 2
#endif

#ifndef THREAD_BINDING
#define THREAD_BINDING spread
#endif

std::vector<Track> simtracks;

std::vector<Track> smat_tracks;
std::vector<Track> plex_tracks;

//==============================================================================

MkFitter *g_mkfp;

const int Nhits = 10; // XXXXX ARGH !!!! What if there's a missing / double layer?

const int Nloop = 100;

//==============================================================================

long64 single_run(int                 n_tracks,
                  MkFitter           *mkfp,
                  std::vector<Track> &trk_in,
                  std::vector<Track> &trk_out)
{
  int theEnd = n_tracks;

  for (int itrack = 0; itrack < theEnd; itrack += NN)
  {
    int end = std::min(itrack + NN, theEnd);

    mkfp->InputTracksAndHits(trk_in, itrack, end);

    for (int x = 0; x < Nloop; ++x)
    {
      if (x != 0)
      {
        mkfp->InputTracksOnly(trk_in, itrack, end);
      }

      mkfp->FitTracks();
    }

#ifndef NO_ROOT
    mkfp->OutputFittedTracks(trk_out, itrack, end);
#endif
  }

  // For propagateLine
  // return long64(Nloop) * n_tracks * (68 + 306) * Nhits;

  return long64(Nloop) * n_tracks * (1200 + 306) * Nhits;
}

//------------------------------------------------------------------------------

long64 single_run_glob(int n_tracks)
{
  return single_run(n_tracks, g_mkfp, simtracks, plex_tracks);
}

//==============================================================================

void test_matriplex()
{
  int Nmin  = 16;
  int Nmax  = 64 * 1024; // 32 * 1024;

  generateTracks(simtracks, Nmax);

  g_mkfp = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Nhits);

  g_mkfp->CheckAlignment();

  Timing t([&](int n_vec)
           {
             return single_run_glob(n_vec);
           });

  t.print_tuple_header();

  // Warm up
  t.time_loop(Nmax, 4);

  for (int n = Nmin; n <= Nmax; n *= 2)
  {
    // Nloop = std::numeric_limits<int>::max() / (374 * Nhits * n);

    // printf("XXX n=%d, nloop=%d\n", n, Nloop);

    t.auto_time_loop(n, 2);
    
    // t.time_loop(n, 1);

    t.print(n);
  }

  _mm_free(g_mkfp);
}

//==============================================================================

void test_vtune()
{
  int Nmax  = 512;

#pragma omp parallel for num_threads(NUM_THREADS)
  for (int i = 0; i < NUM_THREADS; ++i)
  {
    std::vector<Track> sim_trk;
    std::vector<Track> rec_trk;

    generateTracks(sim_trk, Nmax);
    rec_trk.resize(Nmax);

    MkFitter *mf = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Nhits);

    mf->CheckAlignment();

    for (int i = 0 ; i < 200; ++i) 
    {
      single_run(Nmax, mf, sim_trk, rec_trk);
    }

    _mm_free(mf);
  }
}

//==============================================================================

void test_standard()
{
  int  Ntracks  = 64 * 1024;
  bool saveTree = false;

  generateTracks(simtracks, Ntracks);

  double tmp, tsm;

  smat_tracks.reserve(simtracks.size());
  tsm = runFittingTest(simtracks, smat_tracks);

  plex_tracks.resize(simtracks.size());
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
  // test_matriplex();

  test_vtune();

  // test_standard();

  return 0;
}
