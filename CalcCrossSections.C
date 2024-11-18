#include "TFile.h"
#include "TTree.h"

// Make it easier to select components of the arrays (x,y,z,t components of 4 momenta and magnitude)
enum comp {z,x,y,t,m};

// Normalise each vertical strip of a 2D histogram 
// Inputs: p_spline gives the total cross section for C nuclei, p_total is the total number of events
// on C nuclei in the ntuple for each neutrino energy, p_hist is the distribution of events in the ntuple in
// neutrino energy and some other variable on the y axis
void NormaliseToSpline(TGraph* p_spline,TH1D* p_total,TH2D* p_hist){
  for(int i_e=1;i_e<p_hist->GetNbinsX()+1;i_e++){
    double xsec = p_spline->Eval(p_hist->GetXaxis()->GetBinCenter(i_e)); 
    double total = p_total->GetBinContent(i_e);
    if(total == 0) continue;
    for(int i_v=1;i_v<p_hist->GetNbinsY()+1;i_v++){
      p_hist->SetBinContent(i_e,i_v,xsec*p_hist->GetBinContent(i_e,i_v)/total);
      //std::cout << p_hist->GetBinContent(i_e,i_v) << std::endl;
    }
  }
}

void CalcCrossSections(){

  // Load the total cross section spline - only use the bound nucleons for now
  TFile* p_fxsec = TFile::Open("NuanceSplines.root");
  TGraph* p_num_c = static_cast<TGraph*>(p_fxsec->Get("num_c"));
  gROOT->cd();
  p_fxsec->Close();

  TFile* p_fin = TFile::Open("nuance_v3_may07_20070507_numu_ma1.35_kappa1.007_all.root");
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
  //Float_t         x;
  //Float_t         y;
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
  //p_tin->SetBranchAddress("x", &x);
  //p_tin->SetBranchAddress("y", &y);
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
  TH1D* p_nevents = new TH1D("nevents",";Neutrino Energy (GeV);N Events",40,0.0,2.0);

  // Setup histograms - to normalise using the xsec spline we always set the x axis to be 
  // in bins of neutrino energy and the y axis as whatever other varuable we want
  TH2D* p_leptonmomentum_cc = new TH2D("muonmomentum_cc","CC Inclusive;Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
  TH2D* p_leptonmomentum_nc = new TH2D("muonmomentum_nc","NC Inclusive;Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
  TH2D* p_leptoncostheta_cc = new TH2D("muoncostheta_cc","CC Inclusive;Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);
  TH2D* p_leptoncostheta_nc = new TH2D("muoncostheta_nc","NC Inclusive;Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);

  std::map<int,TH2D*> m_ch_leptonmomentum;
  std::map<int,TH2D*> m_ch_leptoncostheta;
  
  for(Long64_t ievent=0;ievent<c_nevents;ievent++){
    p_tin->GetEntry(ievent);

    if(!bound) continue; // only look at interactions on C nuclei for now

    p_nevents->Fill(pneutrino[t]/1e3);

    if(cc) p_leptonmomentum_cc->Fill(pneutrino[t]/1e3,plepton[0][m]/1e3);
    else p_leptonmomentum_nc->Fill(pneutrino[t]/1e3,plepton[0][m]/1e3);

    double costheta = plepton[0][z]/plepton[0][m];
    //std::cout <<  plepton[0][z] << "  " << plepton[0][m] << "  " <<  costheta << std::endl;
    if(cc) p_leptoncostheta_cc->Fill(pneutrino[t]/1e3,costheta);
    else p_leptoncostheta_nc->Fill(pneutrino[t]/1e3,costheta);

    if(m_ch_leptonmomentum.find(channel) == m_ch_leptonmomentum.end()){
      m_ch_leptonmomentum[channel] = new TH2D(Form("muonmomentum_%i",channel),";Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
      m_ch_leptoncostheta[channel] = new TH2D(Form("muoncostheta_%i",channel),";Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);
    }

    m_ch_leptonmomentum[channel]->Fill(pneutrino[t]/1e3,plepton[0][m]/1e3);
    m_ch_leptoncostheta[channel]->Fill(pneutrino[t]/1e3,costheta);

  }

  NormaliseToSpline(p_num_c,p_nevents,p_leptonmomentum_cc);
  NormaliseToSpline(p_num_c,p_nevents,p_leptonmomentum_nc);
  NormaliseToSpline(p_num_c,p_nevents,p_leptoncostheta_cc);
  NormaliseToSpline(p_num_c,p_nevents,p_leptoncostheta_nc);

  TFile* p_fout = new TFile("NuanceCrossSections.root","RECREATE");
  p_fout->cd();
  
  for(std::map<int,TH2D*>::iterator it = m_ch_leptonmomentum.begin();it != m_ch_leptonmomentum.end();it++){
    NormaliseToSpline(p_num_c,p_nevents,it->second);
    it->second->Write(Form("leptonmomentum_%i",it->first));
  }
  for(std::map<int,TH2D*>::iterator it = m_ch_leptoncostheta.begin();it != m_ch_leptoncostheta.end();it++){
    NormaliseToSpline(p_num_c,p_nevents,it->second);
    it->second->Write(Form("leptoncostheta_%i",it->first));
  }

  p_leptonmomentum_cc->Write("leptonmomentum_cc");
  p_leptonmomentum_nc->Write("leptonmomentum_nc");
  p_leptoncostheta_cc->Write("leptoncostheta_cc");
  p_leptoncostheta_nc->Write("leptoncostheta_nc");

  TCanvas* p_canvas = new TCanvas("c","c");

  p_leptonmomentum_cc->Draw("colz");
  p_leptonmomentum_cc->SetStats(0);
  p_canvas->Print("NUANCE_leptonmomentum_cc.png");
  p_canvas->Clear();  

  p_leptonmomentum_nc->Draw("colz");
  p_leptonmomentum_nc->SetStats(0);
  p_canvas->Print("NUANCE_leptonmomentum_nc.png");
  p_canvas->Clear();  
  
  p_leptoncostheta_cc->Draw("colz");
  p_leptoncostheta_cc->SetStats(0);
  p_canvas->Print("NUANCE_leptoncostheta_cc.png");
  p_canvas->Clear();  

  p_leptoncostheta_nc->Draw("colz");
  p_leptoncostheta_nc->SetStats(0);
  p_canvas->Print("NUANCE_leptoncostheta_nc.png");
  p_canvas->Clear();  
  
  p_fout->Close();
  p_fin->Close();
}
