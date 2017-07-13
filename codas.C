#include "Plotter.C"

void codas()
{
  printf("CoDaS plotting macros loaded:\n"
         " - plot_min()\n"
         " - plot_sig()\n"
         " # Following ones require running of additional tests\n"
         " - plot_trig()\n"
         " - plot_vecopt()\n"
         );

}

void plot_min(bool save_p=false)
{
  Plotter p({"sum2",  "sum3", "mul2", "mul3", "div2",  "div3"});
  
  p.XLabel = "N_{array}";
  // p.PlotOpT = p.PlotFlops = false;

  p.do_all("arr", "O3", save_p);
}

void plot_sig(bool save_p=false)
{
  Plotter p({"sum2_sqr",  "sum2_cube", "sum2_quad", "sum2_quint",
             "sum3_sqr",  "sum3_cube"});
  
  p.XLabel = "N_{array}";
  // p.PlotOpT = p.PlotFlops = false;

  p.do_all("arr", "O3", save_p);
}

void plot_trig(bool save_p=false)
{
  Plotter p({"sin2",  "sincos2", "sincos2_tyl4", "sincos2_tyl6", "atan2"});
  
  p.XLabel = "N_{array}";
  // p.PlotOpT = p.PlotFlops = false;

  p.do_all("arr", "O3", save_p);
}

//================================================================

void plot_vecopt(bool save_p=false)
{
  Plotter p({"sum2_cube_O3","sum2_cube_sse4.2", "sum2_cube_avx"});

  p.XLabel = "N_{array}";
  // p.PlotOpT = p.PlotFlops = false;

  p.do_all("arr", "", save_p);
}

//================================================================
// Matriplex plots
//================================================================

void plot_mplex_3_8(bool save_p=false)
{
  Plotter p({"general_3_8", "std_3_8", "intr_3_8",
             "sym_general_3_8", "sym_std_3_8", "sym_intr_3_8"});

  p.XLabel = "N_{array}";
  // p.PlotOpT = p.PlotFlops = false;

  p.do_all("mplx_mult2", "", save_p);
}

void plot_mplex_6_8(bool save_p=false)
{
  Plotter p({"general_6_8", "std_6_8", "intr_6_8",
             "sym_general_6_8", "sym_std_6_8", "sym_intr_6_8"});

  p.XLabel = "N_{array}";
  // p.PlotOpT = p.PlotFlops = false;

  p.do_all("mplx_mult2", "", save_p);
}
