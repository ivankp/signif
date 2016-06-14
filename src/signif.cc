#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <memory>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

#include "binner.hh"
#include "timed_counter.hh"

using namespace std;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

#define h_(name,nbins,lo,hi) TH1D *h_##name = new TH1D(#name,"",nbins,lo,hi);

constexpr pair<double,double> mass_range {105000.,160000.};
constexpr pair<double,double> mass_window{121000.,129000.};

constexpr double len(pair<double,double> p) {
  return p.second - p.first;
}

constexpr double factor = len(mass_window)/(len(mass_range)-len(mass_window));

struct bkg_sig {
  double bkg, sig;
  bkg_sig(): bkg(0), sig(0) { }
  const bkg_sig& operator()(double m, double w) noexcept {
    if (m < mass_window.first || m > mass_window.second) bkg += w;
    else sig += w;
    return *this;
  }
};

struct var {
  string name;
  binner<bkg_sig> bins;
  union { Float_t f; Int_t i; } _x;
  Float_t x() const noexcept {
    if (name[0]=='N') return _x.i;
    else return _x.f;
  }

  var(const string& str) {
    vector<string> tok;
    bool space = true;
    for (char c : str) {
      if (c==' ' || c==',') {
        space = true;
      } else if (space) {
        tok.emplace_back(1,c);
        space = false;
      } else tok.back() += c;
    }
    name = tok.front();
    vector<double> edges(tok.size()-1);
    for (unsigned i=1; i<tok.size(); ++i)
      edges[i-1] = std::stod(tok[i]);
    bins.init(edges.begin(),edges.end());
  }
};

int main(int argc, char* argv[])
{
  const double lumi = 6.;

  vector<unique_ptr<var>> vars;

  ifstream binsf("bins.txt");
  if (binsf.is_open()) {
    string line;
    while ( getline(binsf,line) ) {
      // cout << line << '\n';
      if (line.size()==0 || line[0]=='#') continue;
      vars.emplace_back(new var(line));
    }
    binsf.close();
  } else {
    cout << "Unable to open bins.txt" << endl;
    return 1;
  }

  Char_t isPassed;
  Float_t eff, weight, m_yy;

  timed_counter counter;
  for (int i=1; i<argc; ++i) {
    TFile* file = new TFile(argv[i],"read");
    if (file->IsZombie()) return 1;
    cout << file->GetName() << endl;

    TTree* tree = (TTree*)file->Get("CollectionTree");

    tree->SetBranchAddress("HGamEventInfoAuxDyn.crossSectionBRfilterEff", &eff);
    tree->SetBranchAddress("HGamEventInfoAuxDyn.weight", &weight);
    tree->SetBranchAddress("HGamEventInfoAuxDyn.isPassed", &isPassed);

    tree->SetBranchAddress("HGamEventInfoAuxDyn.m_yy", &m_yy);

    for (auto& v : vars)
      tree->SetBranchAddress(
        ("HGamEventInfoAuxDyn."+v->name).c_str(),
        reinterpret_cast<void*>(&(v->_x))
      );

    for (Long64_t ent=0, n=tree->GetEntries(); ent<n; ++ent) {
      counter(ent);
      tree->GetEntry(ent);

      if (!isPassed) continue;

      if (m_yy < mass_range.first || m_yy > mass_range.second)
        continue;

      const double w = eff*weight*lumi;

      for (auto& v : vars)
        v->bins.fill(v->x(),m_yy,w);
    }
    cout << endl;

    file->Close();
    delete file;
  }

  TFile* file = new TFile("sig.root","recreate");

  vector<TH1*> hists;
  hists.reserve(vars.size());

  for (const auto& v : vars) {
    TH1 *h;
    if (v->name[0]=='N') {
      std::vector<double> edges(v->bins.edges());
      edges.back() = edges[edges.size()-2]+1;
      h = new TH1D(v->name.c_str(),"",v->bins.nbins(),edges.data());
    } else if (v->name[0]=='p' && v->name[1]=='T') {
      std::vector<double> edges(v->bins.edges());
      for (auto& e : edges) e /= 1e3;
      h = new TH1D(v->name.c_str(),"",v->bins.nbins(),edges.data());
    } else {
      h = new TH1D(v->name.c_str(),"",v->bins.nbins(),v->bins.edges().data());
    }
    hists.push_back(h);

    cout << v->name << endl;
    const unsigned n = v->bins.nbins();
    const unsigned w = log10(v->bins.edges().back())+1;
    for (unsigned i=1; i<=n; ++i) {
      cout <<'['<<setw(w)<< v->bins.ledge(i)
           <<','<<setw(w)<< v->bins.redge(i) <<"): "
           << v->bins[i].sig << "  "
           << v->bins[i].bkg << "  ";
      // double sig = v->bins[i].sig / sqrt(
      //   v->bins[i].sig + factor * v->bins[i].bkg );
      double sig = v->bins[i].sig / sqrt( factor * v->bins[i].bkg );
      cout << sig << endl;

      h->SetBinContent(i,sig);
    }
    cout << endl;
  }

  file->Write();
  file->Close();
  delete file;

  return 0;
}
