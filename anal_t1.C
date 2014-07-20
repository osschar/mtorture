// Run this in compiled mode ... the string handling got too complicated :(

#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"

//==============================================================================

const int MaxNN = 6;
      int NN    = 3;

// Vec / no vec comparison; t1-0:
// TString  Gsuff[] = { "host_avx", "mic_mic", "host_novec", "mic_novec", "host_nosimd", "mic_nosimd" };

// Host & mic; t1-1:
// TString  Gsuff[] = { "host_O3", "mic_O3", "cube_host_O3", "cube_mic_O3",  "quint_host_O3", "quint_mic_O3" };

// O3 vs. O2; t1-1:
// TString  Gsuff[] = { "host_O3", "host_O2", "cube_host_O3", "cube_host_O2",  "quint_host_O3", "quint_host_O2" };
// TString  Gsuff[] = { "mic_O3", "mic_O2", "cube_mic_O3", "cube_mic_O2",  "quint_mic_O3", "quint_mic_O2" };

// sum2 (with N=5) and sum3 (with N=3); t1-1:
// TString  Gsuff[] = { "", "sqr", "cube", "quad", "quint" };

// all for SameArray test; t1-2:
// TString  Gsuff[] = {"sum2_cube_sa",  "sum2_cube",
//                     "sum2_quint_sa", "sum2_quint",
//                     "sum3_cube_sa",  "sum3_cube"    };

// t2-0:
// TString  Gsuff[] = {"sum2_cube_host",  "sum2_cube_mic",
//                     "sum2_quint_host", "sum2_quint_mic",
//                     "sum3_cube_host",  "sum3_cube_mic"  };

// t3-0:
// do_all("mult2", "host", save_p); do_all("mult2", "mic", save_p);
TString  Gsuff[] = { "3_8", "3_16", "3_32", "6_8", "6_16", "6_32" };


//------------------------------------------------------------------------------

TString  XLabel = "N_{array}";

bool     PlotVUt = true, PlotOpT = true, PlotFlops = true;

int      Gmsty[] = { 2, 5, 2, 5, 2, 5};
int      Gmcol[] = { kRed+2, kBlue+2, kOrange+4, kCyan+2, kYellow+2, kGreen+2 };
int      Glcol[] = { kRed, kBlue, kOrange, kCyan, kYellow, kGreen };

TTree   *T[MaxNN]   = { 0 };

//==============================================================================

TString join2(const TString& a, const TString& b, const TString& sep="_")
{
  TString j = a;
  if ( ! a.IsNull() && ! b.IsNull())
    j += sep;
  j += b;

  return j;
}

TString join3(const TString& a, const TString& b, const TString& c, const TString& sep="_")
{
  return join2(join2(a, b, sep), c, sep);
}

void load_trees(const TString& pref, const TString& post)
{
  for (int i = 0; i < NN; ++i)
  {
    TString name = join3(pref, Gsuff[i], post) + ".rt";

    T[i] = new TTree;
    T[i]->ReadFile(name);
  }
}

//==============================================================================

TMultiGraph* plot_graph(const TString& varexp,
                        const TString& name, const TString& title,
                        bool save_p=false)
{
  TCanvas *c = new TCanvas(name, name);

  TMultiGraph *mg = new TMultiGraph;
  mg->SetTitle(title);

  TGraph *g;
  int n;

  for (int i = 0; i < NN; ++i)
  {
    n = T[i]->Draw(varexp, "", "goff");
    g = new TGraph(n, T[i]->GetV1(), T[i]->GetV2());

    TString nt(Gsuff[i]); if (nt.IsNull()) nt = " ";
    g->SetName(nt);
    g->SetTitle(nt);

    g->SetMarkerStyle(Gmsty[i]);
    g->SetMarkerColor(Gmcol[i]);
    g->SetLineColor(Glcol[i]);
    g->SetFillColor(kWhite);
    mg->Add(g);
  }

  gPad->SetLogx(true);

  mg->Draw("apl");

  gPad->BuildLegend(0.8, 0.8, 0.99, 0.99);

  if (save_p)
  {
    c->SaveAs(TString::Format("%s.png", name.Data()));
  }

  return mg;
}

//==============================================================================

void plot_all_graphs(const TString& name, const TString& post, bool save_p=false)
{
  TString pref = join2(name, post);

  if (PlotVUt)
    plot_graph("NVec:VecUt", pref + "_VecUt" ,
               join3("Vector utilization", XLabel, "Utilization", ";"), save_p);

  if (PlotOpT)
    plot_graph("NVec:OpT", pref + "_OpT",
               join3("Operations per tick", XLabel, "N_{ops}", ";"), save_p);

  if (PlotFlops)
    plot_graph("NVec:Gflops", pref + "_Gflops",
               join3("Gflops", XLabel, "Gflops", ";"), save_p);
}

//------------------------------------------------------------------------------

void do_all(const TString& name, const TString& post, bool save_p=false)
{
  load_trees(name, post);
  plot_all_graphs(name, post, save_p);
}

//------------------------------------------------------------------------------

void anal_t1(bool save_p=false)
{
  XLabel = "Matriplex vector size";
  PlotOpT = PlotFlops = false;

  NN = 6;
  do_all("mult2_sym", "host", save_p);
  do_all("mult2_sym", "mic", save_p);

  NN = 3;
  do_all("inv_cramer_sym", "host", save_p);
  do_all("inv_cramer_sym", "mic", save_p);
  do_all("inv_cholesky_sym", "host", save_p);
  do_all("inv_cholesky_sym", "mic", save_p);
}
