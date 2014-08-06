#include "MatriplexCommon.h"

#include "fittest.h"

int main()
{
  int  Ntracks  = 5000;
  bool saveTree = false;

  std::vector<Track> simtracks;
  generateTracks(simtracks, Ntracks);

  std::vector<Track> smat_tracks; smat_tracks.reserve(simtracks.size());
  std::vector<Track> plex_tracks; plex_tracks.resize (simtracks.size());


  double tmp, tsm;

  tsm = runFittingTest(simtracks, smat_tracks);

  tmp = runFittingTestPlex(simtracks, plex_tracks);

  printf("SMatrix = %.3f   Matriplex = %.3f   ---   SM/MP = %.3f\n", tsm, tmp, tsm / tmp);

#ifndef NO_ROOT
  make_validation_tree("validation-smat.root", simtracks, smat_tracks);
  make_validation_tree("validation-plex.root", simtracks, plex_tracks);
#endif

  return 0;
}
