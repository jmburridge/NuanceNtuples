#include "TFile.h"
#include "TTree.h"

// Make it easier to select components of the arrays (x,y,z,t components of 4 momenta and magnitude)
enum comp {z,x,y,t,m};

// Normalise each vertical strip of a 2D histogram 
// Inputs: p_spline gives the total cross section for C nuclei, p_total is the total number of events
// on C nuclei in the ntuple for each neutrino energy, p_hist is the distribution of events in the ntuple in
// neutrino energy and some other variable on the y axis
// modified it to an ambiguous 'target' particle: tgt. this can be either a bound  proton or a neutron.
void NormaliseToSpline(TGraph* p_spline_tgt, TH1D* p_total_tgt, TH2D* p_hist_tgt){
  for(int i_e=1;i_e<p_hist_tgt->GetNbinsX()+1;i_e++){
     double xsec_tgt = p_spline_tgt->Eval(p_hist_tgt->GetXaxis()->GetBinCenter(i_e)); 
     double total_tgt = p_total_tgt->GetBinContent(i_e);
    if(total_tgt == 0) continue;
    for(int i_v=1;i_v<p_hist_tgt->GetNbinsY()+1;i_v++){
      p_hist_tgt->SetBinContent(i_e,i_v,xsec_tgt*p_hist_tgt->GetBinContent(i_e,i_v)/total_tgt);
 //std::cout << p_hist->GetBinContent(i_e,i_v) << std::endl;
    }
  }
}

void CalcCrossSections(){

  std::string rootdir = "/exp/uboone/data/users/jburridg/Nuance/";

  // Load the total cross section spline - only use the bound nucleons for now
  TFile* p_fxsec = TFile::Open("NuanceSplines.root");
//  TGraph* p_num_c = static_cast<TGraph*>(p_fxsec->Get("num_c"));
  //Adding in the bound proton and neutron xsec splines too.
  //This will allow us to separate the C nuclei into protons and neutrons. 
  TGraph* p_num_nbcc = static_cast<TGraph*>(p_fxsec->Get("num_nbcc"));
  TGraph* p_num_nbnc = static_cast<TGraph*>(p_fxsec->Get("num_nbnc"));
  TGraph* p_num_pbcc = static_cast<TGraph*>(p_fxsec->Get("num_pbcc"));
  TGraph* p_num_pbnc = static_cast<TGraph*>(p_fxsec->Get("num_pbnc"));
  // Ensures the graph object is accessible after closing the file. 
  gROOT->cd();
  // closes file.
  p_fxsec->Close();
  
  //Open the file and find the tree
  TFile* p_fin = TFile::Open((rootdir + "NUANCE/" + "nuance_v3_may07_20070507_numu_ma1.35_kappa1.007_all.root").c_str());
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
  TH1D* p_nevents_p = new TH1D("nevents_p",";Neutrino Energy (GeV);N Events",40,0.0,2.0);
  TH1D* p_nevents_n = new TH1D("nevents_n",";Neutrino Energy (GeV);N Events",40,0.0,2.0);
  // Setup histograms - to normalise using the xsec spline we always set the x axis to be 
  // in bins of neutrino energy and the y axis as whatever other variable we want
  // create more histograms here for protons and neutrons. 

  //define histograms for lepton momentum   
  TH2D* p_leptonmomentum_cc_p = new TH2D("muonmomentum_cc","CC Proton;Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
  TH2D* p_leptonmomentum_cc_n = new TH2D("muonmomentum_cc","CC Neutron Target;Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);

  TH2D* p_leptonmomentum_nc_p = new TH2D("muonmomentum_nc","NC Proton Target;Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
  TH2D* p_leptonmomentum_nc_n = new TH2D("muonmomentum_nc","NC Neutron Target;Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
 
  // Define Histograms for lepton costheta
  TH2D* p_leptoncostheta_cc_p = new TH2D("muoncostheta_cc","CC Proton Target;Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);
 TH2D* p_leptoncostheta_cc_n = new TH2D("muoncostheta_cc","CC Neutron Target;Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);

  TH2D* p_leptoncostheta_nc_p = new TH2D("muoncostheta_nc","NC Proton Target;Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);
  TH2D* p_leptoncostheta_nc_n = new TH2D("muoncostheta_nc","NC Neutron Target;Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);

  
  std::map<std::pair<int, int>, TH2D*> m_ch_leptonmomentum;  // Using pair to hold channel and target type
  std::map<std::pair<int, int>, TH2D*> m_ch_leptoncostheta;

  for(Long64_t ievent=0;ievent<c_nevents;ievent++){
    p_tin->GetEntry(ievent);

    if(!bound) continue; // only look at interactions on C nuclei for now 
    p_nevents->Fill(pneutrino[t]/1e3);


    int nucleon_type = (target == 2212) ? 2212 : 2112;  // Distinguish proton (2212) and neutron (2112)
    std::pair<int, int> key = std::make_pair(channel, nucleon_type);

    double costheta = plepton[0][z]/plepton[0][m];

    //adding boolean values for proton and neutron traget identification
    bool target_is_proton = (target = 2122);
    bool target_is_neutron = (target = 2112); //PDG ID codes 

    // adding nested conditional statements to distinguish between proton and neutron targets. 
    if (target_is_proton){
        p_nevents_p->Fill(pneutrino[t]/1e3);
        if(cc) {
	       
		p_leptonmomentum_cc_p->Fill(pneutrino[t]/1e3,plepton[0][m]/1e3);
                p_leptoncostheta_cc_p->Fill(pneutrino[t]/1e3,costheta);
	}
	else {       
		p_leptonmomentum_nc_p->Fill(pneutrino[t]/1e3,plepton[0][m]/1e3);
                p_leptoncostheta_nc_p->Fill(pneutrino[t]/1e3,costheta);
	}
    }
    else if (target_is_neutron) {
       p_nevents_n->Fill(pneutrino[t]/1e3);
       if(cc) {
	      p_leptonmomentum_cc_n->Fill(pneutrino[t]/1e3,plepton[0][m]/1e3);
              p_leptoncostheta_cc_n->Fill(pneutrino[t]/1e3,costheta);
       }
       else {
	       p_leptonmomentum_nc_n->Fill(pneutrino[t]/1e3,plepton[0][m]/1e3);
               p_leptoncostheta_nc_n->Fill(pneutrino[t]/1e3,costheta);

       }
    } 
    //std::cout <<  plepton[0][z] << "  " << plepton[0][m] << "  " <<  costheta << std::endl;


     if (m_ch_leptonmomentum.find(key) == m_ch_leptonmomentum.end()) {
        // Create new histograms for each channel and nucleon type
        m_ch_leptonmomentum[key] = new TH2D(Form("leptonmomentum_%i_%i", channel, nucleon_type), ";Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)", 40, 0.0, 2.0, 40, 0.0, 2.0);
        m_ch_leptoncostheta[key] = new TH2D(Form("leptoncostheta_%i_%i", channel, nucleon_type), ";Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)", 40, 0.0, 2.0, 40, -1.0, 1.0);
    }

    m_ch_leptonmomentum[key]->Fill(pneutrino[t]/1e3,plepton[0][m]/1e3);
    m_ch_leptoncostheta[key]->Fill(pneutrino[t]/1e3,costheta);

  }

  NormaliseToSpline(p_num_pbcc,p_nevents_p,p_leptonmomentum_cc_p);
  NormaliseToSpline(p_num_pbnc,p_nevents_p,p_leptonmomentum_nc_p); //p_nevents_p needs to be seprated into nc and cc? 
  NormaliseToSpline(p_num_nbcc,p_nevents_n,p_leptonmomentum_cc_n);
  NormaliseToSpline(p_num_nbnc,p_nevents_n,p_leptonmomentum_nc_n);
  NormaliseToSpline(p_num_pbcc,p_nevents_p,p_leptoncostheta_cc_p);
  NormaliseToSpline(p_num_pbnc,p_nevents_p,p_leptoncostheta_nc_p);
  NormaliseToSpline(p_num_nbcc,p_nevents_n,p_leptoncostheta_cc_n);
  NormaliseToSpline(p_num_nbnc,p_nevents_n,p_leptoncostheta_nc_n);

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

