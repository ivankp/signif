#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <utility>

#include <TFile.h>
#include <TTree.h>
// #include <TH1.h>

#include "binner.hh"
#include "timed_counter.hh"

using namespace std;

#define h_(name,nbins,lo,hi) TH1D *h_##name = new TH1D(#name,"",nbins,lo,hi);

const pair<double,double> mass_window{121.,129.};

struct bkg_sig {
  double bkg, sig;
  bkg_sig(): bkg(0), sig(0) { }
  const bkg_sig& operator()(double m, double w) noexcept {
    if (m < mass_window.first || m > mass_window.second) bkg += w;
    else sig += w;
    return *this;
  }
};

int main(int argc, char* argv[])
{
  const double lumi = 6.;

  // TFile* out = new TFile("out.root","recreate");

  // h_(m_yy,100,105,160)

  binner<bkg_sig> pT_yy_bs({0,20,100,200});

  Float_t eff, weight, m_yy, pT_yy;

  timed_counter counter;
  for (int i=1; i<argc; ++i) {
    TFile* file = new TFile(argv[i],"read");
    if (file->IsZombie()) return 1;
    cout << file->GetName() << endl;

    TTree* tree = (TTree*)file->Get("CollectionTree");

    tree->SetBranchAddress("HGamEventInfoAuxDyn.crossSectionBRfilterEff", &eff);
    tree->SetBranchAddress("HGamEventInfoAuxDyn.weight", &weight);

    tree->SetBranchAddress("HGamEventInfoAuxDyn.m_yy", &m_yy);
    tree->SetBranchAddress("HGamEventInfoAuxDyn.pT_yy", &pT_yy);

    for (Long64_t ent=0, n=100000/*tree->GetEntries()*/; ent<n; ++ent) {
      counter(ent);
      tree->GetEntry(ent);
      m_yy /= 1e3;

      const double w = eff*weight*lumi;

      // h_m_yy->Fill(m_yy,w);
      pT_yy_bs.fill(pT_yy/1e3,m_yy,w);
    }
    cout << endl;

    file->Close();
    delete file;
  }

  for (unsigned i=1, n=pT_yy_bs.nbins(); i<=n; ++i)
    cout <<'['<<setw(3)<< pT_yy_bs.ledge(i)
         <<','<<setw(3)<< pT_yy_bs.redge(i) <<"): "
         << pT_yy_bs[i].sig << '\t'
         << pT_yy_bs[i].bkg << endl;

  // out->Write();
  // out->Close();
  // delete out;

  return 0;
}
