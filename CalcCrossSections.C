#include "TFile.h"
#include "TTree.h"

// estimated number of C atoms in miniboone
// TODO: Find out the fraction of hydrogen/other elements and the exact mass of the
// simulation volume
const double c_ntargs = 800*1e3/0.012*6.20e23;

// Make it easier to select components of the arrays (x,y,z,t components of 4 momenta and magnitude)
enum comp {x,y,z,t,m};

// Normalise each vertical strip of a 2D histogram 
void NormaliseToSpline(TGraph* p_spline,TH1D* p_total,TH2D* p_hist){

  for(int i_e=1;i_e<p_hist->GetNbinsX()+1;i_e++){
    double xsec = p_spline->Eval(p_hist->GetXaxis()->GetBinCenter(i_e)); 
    double total = p_total->GetBinContent(i_e);

    if(total == 0) continue;
    //double integral = p_hist->Integral(i_e,i_e,1,p_hist->GetNbinsY()); 
    //std::cout << "ratio=" <<  integral/total << "=" << integral << "/" << total << std::endl;

    for(int i_v=1;i_v<p_hist->GetNbinsY();i_v++){
      p_hist->SetBinContent(i_e,i_v,xsec*p_hist->GetBinContent(i_e,i_v)/total);
    }
    
  }

}

void CalcCrossSections(){

  // Load the total cross section spline - only use the bound nucleons for now
  TFile* p_fxsec = TFile::Open("NuanceSplines.root");
  TGraph* p_num_c = static_cast<TGraph*>(p_fxsec->Get("num_c"));
  gROOT->cd();
  p_fxsec->Close();

  TFile* p_fin = TFile::Open("nuance_v3_may07_20070507_numu_ma1.35_kappa1.007_1.root");
  TTree* p_tin = static_cast<TTree*>(p_fin->Get("h3"));

  // Setup branches

  UChar_t         cc;
  UChar_t         bound;
  Int_t           event;
  Int_t           neutrino;
  Int_t           target;
  Int_t           iniq;
  Int_t           finq;
  Int_t           lepton0;
  Float_t         polar;
  Int_t           channel;
  Float_t         qsq;
  Float_t         w;
  Float_t         x;
  Float_t         y;
  Float_t         pneutrino[4];
  Float_t         ptarg[5];
  Float_t         vertex[4];
  Float_t         start[4];
  Float_t         depth;
  Float_t         flux;
  Int_t           n_leptons;
  Float_t         pltot[5];
  Int_t           lepton[10];   //[n_leptons]
  Float_t         plepton[10][5];   //[n_leptons]
  Int_t           n_hadrons;
  Float_t         phtot[5];
  Int_t           hadron[50];   //[n_hadrons]
  Float_t         phadron[50][5];   //[n_hadrons]

  p_tin->SetBranchAddress("cc", &cc);
  p_tin->SetBranchAddress("bound", &bound);
  p_tin->SetBranchAddress("event", &event);
  p_tin->SetBranchAddress("neutrino", &neutrino);
  p_tin->SetBranchAddress("target", &target);
  p_tin->SetBranchAddress("iniq", &iniq);
  p_tin->SetBranchAddress("finq", &finq);
  p_tin->SetBranchAddress("lepton0", &lepton0);
  p_tin->SetBranchAddress("polar", &polar);
  p_tin->SetBranchAddress("channel", &channel);
  p_tin->SetBranchAddress("qsq", &qsq);
  p_tin->SetBranchAddress("w", &w);
  p_tin->SetBranchAddress("x", &x);
  p_tin->SetBranchAddress("y", &y);
  p_tin->SetBranchAddress("p_neutrino", pneutrino);
  p_tin->SetBranchAddress("p_targ", ptarg);
  p_tin->SetBranchAddress("vertex", vertex);
  p_tin->SetBranchAddress("start", start);
  p_tin->SetBranchAddress("depth", &depth);
  p_tin->SetBranchAddress("flux", &flux);
  p_tin->SetBranchAddress("n_leptons", &n_leptons);
  p_tin->SetBranchAddress("p_ltot", pltot);
  p_tin->SetBranchAddress("lepton", lepton);
  p_tin->SetBranchAddress("p_lepton", plepton);
  p_tin->SetBranchAddress("n_hadrons", &n_hadrons);
  p_tin->SetBranchAddress("p_htot", phtot);
  p_tin->SetBranchAddress("hadron", hadron);
  p_tin->SetBranchAddress("p_hadron", phadron);

  const Long64_t c_nevents = p_tin->GetEntries();
  //std::cout << "Tree has " << c_nevents << " events" << std::endl;

  // Record the total number of events generated for each neutrino energy
  TH1D* p_nevents = new TH1D("nevents",";Neutrino Energy (GeV);N Events",20,0.0,2.0);

  // Setup histograms - to normalise using the xsec spline we always set the x axis to be 
  // in bins of neutrino energy and the y axis as whatever other varuable we want
  TH2D* p_leptonmomentum_cc = new TH2D("muonmomentum_cc","CC Inclusive;Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",20,0.0,2.0,20,0.0,2.0);
  TH2D* p_leptonmomentum_nc = new TH2D("muonmomentum_nc","NC Inclusive;Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",20,0.0,2.0,20,0.0,2.0);

  
  for(Long64_t ievent=0;ievent<c_nevents;ievent++){
    p_tin->GetEntry(ievent);

    if(!bound) continue; // only look at interactions on C nuclei for now

    p_nevents->Fill(pneutrino[t]/1e3);

    if(cc) p_leptonmomentum_cc->Fill(pneutrino[t]/1e3,plepton[0][m]/1e3);
    else p_leptonmomentum_nc->Fill(pneutrino[t]/1e3,plepton[0][m]/1e3);

  }

  NormaliseToSpline(p_num_c,p_nevents,p_leptonmomentum_cc);
  NormaliseToSpline(p_num_c,p_nevents,p_leptonmomentum_nc);
  
  TCanvas* pcanvas = new TCanvas("c","c");
  pcanvas->Print("test.png");
  

}
