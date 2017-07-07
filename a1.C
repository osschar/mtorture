#include "Plotter.C"

// Vec / no vec comparison; t1-0:
// { "host_avx", "mic_mic", "host_novec", "mic_novec", "host_nosimd", "mic_nosimd" };

// Host & mic; t1-1:
// { "host_O3", "mic_O3", "cube_host_O3", "cube_mic_O3",  "quint_host_O3", "quint_mic_O3" };


// Host & mic; t1-1:
// { "host_O3", "host_O2", "cube_host_O3", "cube_host_O2",  "quint_host_O3", "quint_host_O2" };
// { "mic_O3", "mic_O2", "cube_mic_O3", "cube_mic_O2",  "quint_mic_O3", "quint_mic_O2" };

// SameArray test; t1-2:
// { "sum2_cube_sa",  "sum2_cube", "sum2_quint_sa", "sum2_quint", "sum3_cube_sa", "sum3_cube" };


void a1(bool save_p=false)
{
  Plotter p({"sum2_sqr",  "sum2_cube", "sum2_quad", "sum2_quint",
             "sum3_sqr",  "sum3_cube"});
  
  p.XLabel = "N_{array}";
  // p.PlotOpT = p.PlotFlops = false;

  p.do_all("arr", "host", save_p);
}
