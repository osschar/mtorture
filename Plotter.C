#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"

//==============================================================================
// Utility
//------------------------------------------------------------------------------

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

//==============================================================================
// Plotter
//------------------------------------------------------------------------------

// Max number of of graphs on one plot. Limited by marker style arrays.

const int MaxNN = 6;

class Plotter
{
public:
  int      NN;
  TTree   *T[MaxNN];

  TString  XLabel = "N_{array}";

  TString  Gsuff[MaxNN];

  int      Gmsty[MaxNN] = { 2, 5, 2, 5, 2, 5};
  int      Gmcol[MaxNN] = { kRed+2, kBlue+2, kOrange+4, kCyan+2, kMagenta+2, kGreen+2 };
  int      Glcol[MaxNN] = { kRed, kBlue, kOrange, kCyan, kMagenta, kGreen };

  float    Gmsize  = 1.5;
  int      Glwidth = 2;

  bool     PlotVUt = true, PlotOpT = true, PlotFlops = true;

  //==============================================================================

  Plotter(int n) : NN(n)
  {
    if (NN > MaxNN) throw std::runtime_error("Plotter only up to six graphs supported ... increseextend color / marker arrays.");
  }

  Plotter(std::initializer_list<const char*> ss) : Plotter(ss.size())
  {
    set_suffs(ss);
  }

  void set_suffs(std::initializer_list<const char*> ss)
  {
    int i = 0;
    for (auto s : ss) Gsuff[i++] = s;
  }

  //==============================================================================

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

      g->SetMarkerSize(Gmsize);
      g->SetLineWidth(Glwidth);
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

  void plot_all_graphs(const TString& pref, const TString& post, bool save_p=false)
  {
    TString outf = join2(pref, post);

    if (PlotVUt)
      plot_graph("NVec:VecUt", outf + "_VecUt" ,
                 join3("Vector utilization", XLabel, "Utilization", ";"), save_p);

    if (PlotOpT)
      plot_graph("NVec:OpT", outf + "_OpT",
                 join3("Operations per tick", XLabel, "N_{ops}", ";"), save_p);

    if (PlotFlops)
      plot_graph("NVec:Gflops", outf + "_Gflops",
                 join3("Gflops", XLabel, "Gflops", ";"), save_p);
  }

  //------------------------------------------------------------------------------

  void do_all(const TString& pref, const TString& post, bool save_p=false)
  {
    load_trees(pref, post);
    plot_all_graphs(pref, post, save_p);
  }

};
