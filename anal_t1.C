#include "TTree.h"
#include "TMultiGraph.h"

//==============================================================================

const int NN = 3;

TTree   *T[NN]   = { 0 };

// Vec / no vec comparison; t1-0:
// TString  Gsuff[] = { "host_avx", "mic_mic", "host_novec", "mic_novec", "host_nosimd", "mic_nosimd" };

// Host & mic; t1-1:
// TString  Gsuff[] = { "host_O3", "mic_O3", "cube_host_O3", "cube_mic_O3",  "quint_host_O3", "quint_mic_O3" };

// O3 vs. O2; t1-1:
// TString  Gsuff[] = { "host_O3", "host_O2", "cube_host_O3", "cube_host_O2",  "quint_host_O3", "quint_host_O2" };
// TString  Gsuff[] = { "mic_O3", "mic_O2", "cube_mic_O3", "cube_mic_O2",  "quint_mic_O3", "quint_mic_O2" };

// sum2 (with N=5) and sum3 (with N=3); t1-1:
TString  Gsuff[] = { "", "sqr", "cube", "quad", "quint" };


// TString  Gsuff[] = { "_O2", "_O3", "cube_O2", "cube_O3" };

int      Gmsty[] = { 2, 5, 2, 5, 2, 5};
int      Gmcol[] = { kRed+2, kBlue+2, kOrange+4, kCyan+2, kYellow+2, kGreen+2 };
int      Glcol[] = { kRed, kBlue, kOrange, kCyan, kYellow, kGreen };

//==============================================================================

TString join2(const TString& a, const TString& b)
{
  TString j = a;
  if ( ! a.IsNull() && ! b.IsNull())
    j += "_";
  j += b;

  return j;
}

TString join3(const TString& a, const TString& b, const TString& c)
{
  return join2(join2(a, b), c);
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
    g->SetName(Gsuff[i]);
    g->SetTitle(Gsuff[i]);
    g->SetMarkerStyle(Gmsty[i]);
    g->SetMarkerColor(Gmcol[i]);
    g->SetLineColor(Glcol[i]);
    g->SetFillColor(kWhite);
    mg->Add(g);
  }

  gPad->SetLogx(true);

  mg->Draw("apl");

  gPad->BuildLegend(0.8, 0.8, 1, 1);

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

  plot_graph("NVec:VecUt", pref + "_VecUt" ,
             "Vector utilization;N_{array};Utilization", save_p);

  plot_graph("NVec:OpT", pref + "_OpT",
             "Operations per tick;N_{array};N_{ops}", save_p);

  plot_graph("NVec:Gflops", pref + "_Gflops",
             "Gflops;N_{array};Gflops", save_p);
}

//------------------------------------------------------------------------------

void do_all(const TString& name, const TString& post, bool save_p=false)
{
  load_trees(name, post);
  plot_all_graphs(name, post, save_p);
}

//------------------------------------------------------------------------------

void anal1_t1(bool save_p=false)
{
  do_all("arr_sum3", "mic_O3", save_p);
  // do_all("arr_sum2_cube", save_p);
}
