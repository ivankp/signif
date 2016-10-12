#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <memory>
#include <cmath>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TKey.h>
#include <TAxis.h>

#include "branches.hh"
#include "binner.hh"
#include "timed_counter.hh"
#include "in.hh"

using namespace std;
namespace po = boost::program_options;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

constexpr pair<double,double> mass_range {105000.,160000.};
constexpr pair<double,double> mass_window{121000.,129000.};

constexpr double len(pair<double,double> p) {
  return p.second - p.first;
}
template <typename T1, typename T2, typename T3>
inline bool in(T1 x, const pair<T2,T3>& p) noexcept {
  return (p.first < x && x < p.second);
}

constexpr double factor = len(mass_window)/(len(mass_range)-len(mass_window));

struct { double in, need, fac; } lumi;

bool mc_file;

class bkg_sig {
  double bkg_tmp, sig_tmp;
public:
  double bkg, sig;
  bkg_sig(): bkg_tmp(0), sig_tmp(0), bkg(0), sig(0) { }
  void operator()(double m, double w) noexcept {
    bool in_window = in(m,mass_window);
    if ( mc_file && in_window ) sig_tmp += w;
    else if ( !mc_file && !in_window ) bkg_tmp += w;
  }
  void merge(double n) {
    if (mc_file) {
      sig += sig_tmp/n;
      sig_tmp = 0.;
    } else {
      bkg += bkg_tmp/n;
      bkg_tmp = 0.;
    }
  }
  double signif() {
    return sig!=0 ? lumi.fac * sig / sqrt(sig + factor * bkg) : 0;
  }
};

struct var {
  string name;
  binner<bkg_sig> bins;
  bool take_abs;
  union { Float_t f; Int_t i; } _x;
  Float_t x() const noexcept {
    Float_t out = (name[0]=='N') ? _x.i : _x.f;
    if (take_abs) out = abs(out);
    return out;
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

    take_abs = in(name,"Dphi_j_j","cosTS_yy","Dy_j_j");
  }

  double signif(unsigned i) {
    return bins[i].signif();
  }

  void merge(double n) {
    for (auto& b : bins.bins()) b.merge(n);
  }
};

int main(int argc, char* argv[])
{
  vector<string> ifname_data, ifname_mc;
  string ofname, cfname, ifname_bins;

  // options ---------------------------------------------------
  try {
    po::options_description desc("Options");
    desc.add_options()
      ("data", po::value(&ifname_data)->multitoken()->required(),
       "input root data files")
      ("mc", po::value(&ifname_mc)->multitoken()->required(),
       "input root Monte Carlo files")
      ("output,o", po::value(&ofname)->required(),
       "output root file")
      ("conf,c", po::value(&cfname),
       "configuration file")
      ("bins,b", po::value(&ifname_bins)->required(),
       "differential variables bins")
      ("lumi.in", po::value(&lumi.in)->default_value(3245.),
       "configuration file")
      ("lumi.need,l", po::value(&lumi.need)->default_value(6000.),
       "configuration file")
    ;

    po::positional_options_description pos;
    pos.add("conf",1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
      .options(desc).positional(pos).run(), vm);
    if (argc == 1) {
      cout << desc << endl;
      return 0;
    }
    if (vm.count("conf")) {
      po::store( po::parse_config_file<char>(
        vm["conf"].as<string>().c_str(), desc), vm);
    }
    po::notify(vm);
  } catch (exception& e) {
    cerr << "\033[31m" << argv[0]
         << " options: " <<  e.what() <<"\033[0m"<< endl;
    return 1;
  }
  // end options ---------------------------------------------------

  vector<pair<const string*,bool>> ifname;
  for (const auto& f : ifname_data) ifname.emplace_back(&f,false);
  for (const auto& f : ifname_mc  ) ifname.emplace_back(&f,true );

  lumi.fac = sqrt(lumi.need/lumi.in);

  bkg_sig inclusive, nj0, njless2;
  vector<unique_ptr<var>> vars;

  ifstream binsf(ifname_bins);
  if (binsf.is_open()) {
    string line;
    while ( getline(binsf,line) ) {
      // cout << line << '\n';
      if (line.size()==0 || line[0]=='#') continue;
      vars.emplace_back(new var(line));
    }
    binsf.close();
  } else {
    cout << "Unable to open " << ifname_bins << endl;
    return 1;
  }

  Char_t isPassed;
  Float_t cs_br_fe, weight, m_yy;
  Int_t *njets = nullptr;

  for (const auto& f : ifname) {
    TFile* file = new TFile(f.first->c_str(),"read");
    if (file->IsZombie()) return 1;

    mc_file = f.second;

    cout << ( mc_file ? "MC:" : "Data:" ) << ' '
         << file->GetName() << endl;

    double n_all = 1;
    if (mc_file) {
      TIter next(file->GetListOfKeys());
      TKey *key;
      while ((key = static_cast<TKey*>(next()))) {
        string name(key->GetName());
        if (name.substr(0,8)!="CutFlow_" ||
            name.substr(name.size()-18)!="_noDalitz_weighted") continue;
        TH1 *h = static_cast<TH1*>(key->ReadObj());
        cout << h->GetName() << endl;
        n_all = h->GetBinContent(3);
        cout << h->GetXaxis()->GetBinLabel(3) << " = " << n_all << endl;
        break;
      }
    }

    TTree* tree = (TTree*)file->Get("CollectionTree");

    branches(tree,
      "HGamEventInfoAuxDyn.weight",   &weight,
      "HGamEventInfoAuxDyn.isPassed", &isPassed,
      "HGamEventInfoAuxDyn.m_yy",     &m_yy
    );

    if (mc_file) branches_set_on(tree,
      "HGamEventInfoAuxDyn.crossSectionBRfilterEff", &cs_br_fe);

    for (auto& v : vars) {
      branches_set_on(tree,
        ("HGamEventInfoAuxDyn."+v->name).c_str(),
        reinterpret_cast<void*>(&(v->_x))
      );
      if (v->name=="N_j_30") {
        njets = &v->_x.i;
      }
    }

    if (!njets) {
      cerr << "N_j_30 branch was not used" << endl;
      return 1;
    }

    for (timed_counter<Long64_t> ent(tree->GetEntries()); ent.ok(); ++ent) {
      tree->GetEntry(ent);

      if (!isPassed) continue;
      if (!in(m_yy,mass_range)) continue;

      if (mc_file) {
        weight *= cs_br_fe*lumi.in;
      }

      inclusive(m_yy,weight);
      for (auto& v : vars) {
        v->bins.fill(v->x(),m_yy,weight);
      }
      if ((*njets) == 0) nj0(m_yy,weight);
      if ((*njets)  < 2) njless2(m_yy,weight);
    }

    inclusive.merge(n_all);
    nj0.merge(n_all);
    njless2.merge(n_all);
    for (auto& v : vars) v->merge(n_all);

    file->Close();
    delete file;
  }

  cout << "============" << endl;
  test(factor)
  cout << "Inclusive" << endl;
  cout << "Signal: " << inclusive.sig << endl;
  cout << "Bkg under signal: " << factor * inclusive.bkg << endl;
  cout << "Bkg in sidebands: " << inclusive.bkg << endl;
  cout << "Significance: " << inclusive.signif() << endl;
  cout << "============" << endl;

  TFile* file = new TFile(ofname.c_str(),"recreate");

  vector<TH1*> hists;
  hists.reserve(vars.size());

  for (const auto& v : vars) {
    std::vector<double>& edges = v->bins.edges();
    if ( std::isinf(edges.back()) ) {
      edges.back() = edges[edges.size()-2] + (edges[1] - edges[0]);
    }
    if ( (v->name[0]=='p' && v->name[1]=='T')
      || (v->name[0]=='m' && v->name[1]=='_') ) {
      for (auto& e : edges) e /= 1e3;
    }

    const unsigned n = v->bins.nbins();
    TH1 *h = new TH1D(v->name.c_str(),"",n,v->bins.edges().data());
    hists.push_back(h);

    cout << v->name << endl;
    unsigned bin = 1;
    if ( v->name.substr(v->name.size()-4)=="_j_j"
      || v->name.substr(v->name.size()-3)=="_jj"
    ) {
      h->SetBinContent(bin++,njless2.signif());
    } else if ( v->name.substr(v->name.size()-3)=="_j1" ) {
      h->SetBinContent(bin++,nj0.signif());
    }

    const unsigned w = log10(v->bins.edges().back())+1;
    for (; bin<=n; ++bin) {
      cout <<'['<<setw(w)<< v->bins.ledge(bin)
           <<','<<setw(w)<< v->bins.redge(bin) <<"): "
           << v->bins[bin].sig << "  "
           << v->bins[bin].bkg << "  ";

      double signif = v->signif(bin);
      cout << signif << endl;
      h->SetBinContent(bin,signif);
    }
    cout << endl;

    TAxis *xa = h->GetXaxis();
    if (v->name[0]=='N') {
      for (unsigned i=1; ; ++i) {
        stringstream ss;
        if (n-i) {
          ss << " = " << i-1;
          xa->SetBinLabel(i,ss.str().c_str());
        } else {
          ss << " #geq " << i-1;
          xa->SetBinLabel(i,ss.str().c_str());
          break;
        }
      }
      xa->SetLabelSize(0.05);
    }
  }

  file->Write();
  file->Close();
  delete file;

  return 0;
}
