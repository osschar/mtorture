#include "Plotter.C"

void codas()
{
  printf("CoDaS plotting macros loaded:\n"
         " - plot_min()\n"
         " - plot_sig()\n"
         " - plot_trig()\n"
         );

}

void plot_min(bool save_p=false)
{
  Plotter p({"sum2",  "sum3", "mul2", "mul3", "div2",  "div3"});
  
  p.XLabel = "N_{array}";
  // p.PlotOpT = p.PlotFlops = false;

  p.do_all("arr", "min", save_p);
}

void plot_sig(bool save_p=false)
{
  Plotter p({"sum2_sqr",  "sum2_cube", "sum2_quad", "sum2_quint",
             "sum3_sqr",  "sum3_cube"});
  
  p.XLabel = "N_{array}";
  // p.PlotOpT = p.PlotFlops = false;

  p.do_all("arr", "sig", save_p);
}

void plot_trig(bool save_p=false)
{
  Plotter p({"sin2",  "sincos2", "sincos2_tyl4", "sincos2_tyl6", "atan2"});
  
  p.XLabel = "N_{array}";
  // p.PlotOpT = p.PlotFlops = false;

  p.do_all("arr", "trig", save_p);
}
