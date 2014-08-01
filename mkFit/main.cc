/*
  g++ -std=c++11 -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
  icc -std=gnu++0x -O3 -openmp -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
*/

#include <iostream>

#include "MatriplexCommon.h"

#include "fittest.h"
//#include "buildtest.h"

int main()
{
  bool saveTree = false;

  double tmp, tsm;

  tsm = runFittingTest(saveTree, 5000);

  tmp = runFittingTestPlex(saveTree, 5000);

  printf("SMatrix = %.3f   Matriplex = %.3f   ---   SM/MP = %.3f\n", tsm, tmp, tsm / tmp);

  // runBuildingTest(saveTree,10);

  return 0;
}
