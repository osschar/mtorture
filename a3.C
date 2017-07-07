#include "Plotter.C"

void a3(bool save_p=false)
{
  Plotter p({ "3_8", "3_16", "3_32", "6_8", "6_16", "6_32" });

  p.XLabel = "Matriplex vector size";
  p.PlotOpT = p.PlotFlops = false;

  do_all("mult2_sym", "host", save_p);
  do_all("mult2_sym", "mic", save_p);

  p.NN = 3; // No inversion for 6x6 matrices

  do_all("inv_cramer_sym", "host", save_p);
  do_all("inv_cramer_sym", "mic", save_p);
  do_all("inv_cholesky_sym", "host", save_p);
  do_all("inv_cholesky_sym", "mic", save_p);
}
