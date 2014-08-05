#ifndef _fittest_
#define _fittest_

#include "Track.h"

void   generateTracks(std::vector<Track>& simtracks, int Ntracks);

double runFittingTest(std::vector<Track>& simtracks, bool saveTree);

#ifndef __APPLE__
double runFittingTestPlex(std::vector<Track>& simtracks, bool saveTree);
#endif

#endif
