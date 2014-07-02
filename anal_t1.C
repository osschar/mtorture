#include "TTree.h"
#include "TMultiGraph.h"

//==============================================================================

const int NN = 4;

TTree   *T[NN]   = { 0 };
TString  Gsuff[] = { "host_avx", "mic_mic", "host_novec", "mic_novec", "host_nosimd", "mic_nosimd" };

int      Gmsty[] = { 2, 5, 2, 5, 2, 5};
int      Gmcol[] = { kRed+2, kBlue+2, kOrange+4, kCyan+2, kYellow+2, kGreen+2 };
int      Glcol[] = { kRed, kBlue, kOrange, kCyan, kYellow, kGreen };

//==============================================================================

void load_trees(const TString& pref)
{
  for (int i = 0; i < NN; ++i)
  {
    TString name = pref + "_" + Gsuff[i] + ".rt";

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

void plot_all_graphs(const TString& name, bool save_p=false)
{
  plot_graph("NVec:VecUt", name + "_VecUt",
             "Vector utilization;N_{array};Utilization", save_p);

  plot_graph("NVec:OpT", name + "_OpT",
             "Operations per tick;N_{array};N_{ops}", save_p);

  plot_graph("NVec:Gflops", name + "_Gflops",
             "Gflops;N_{array};Gflops", save_p);
}

//------------------------------------------------------------------------------

void do_all(const TString& name, bool save_p=false)
{
  load_trees(name);
  plot_all_graphs(name, save_p);
}

//------------------------------------------------------------------------------

void anal1_t1(bool save_p=false)
{
  do_all("arr_sum2", save_p);
  do_all("arr_sum2_cube", save_p);
}
