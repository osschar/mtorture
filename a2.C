#include "Plotter.C"

void a2()
{
  Plotter p({"sum2_cube_host",  "sum2_cube_mic",
             "sum2_quint_host", "sum2_quint_mic",
             "sum3_cube_host",  "sum3_cube_mic"  });

  p.XLabel = "N_{array}";
  // p.PlotOpT = p.PlotFlops = false;

  p.do_all("mt", "", save_p);
}
