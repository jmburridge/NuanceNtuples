#include "TFile.h"
#include "TTree.h"

// Normalise each vertical strip of a 2D histogram 
// Inputs: p_spline gives the total cross section for C nuclei, p_total is the total number of events
// on C nuclei in the ntuple for each neutrino energy, p_hist is the distribution of events in the ntuple in
// neutrino energy and some other variable on the y axis
void NormaliseToSpline(TGraph* p_spline_tgt,TH1D* p_total_tgt,TH2D* p_hist_tgt){
  //std::cout << "Normalising to spline" << std::endl;
  for(int i_e=1;i_e<p_hist_tgt->GetNbinsX()+1;i_e++){
    double xsec = p_spline_tgt->Eval(p_hist_tgt->GetXaxis()->GetBinCenter(i_e)); 
    double total = p_total_tgt->GetBinContent(i_e);
    //std::cout<< total << "," << p_hist_tgt->GetXaxis()->GetBinCenter(i_e)<<", " << xsec << std::endl;
    if(total == 0) continue;
    for(int i_v=1;i_v<p_hist_tgt->GetNbinsY()+1;i_v++){
      std::cout << p_hist_tgt->GetXaxis()->GetBinCenter(i_e) << ", "<<  p_hist_tgt->GetXaxis()->GetBinCenter(i_v) <<", "<< total <<", " <<  p_hist_tgt->GetBinContent(i_e,i_v)<<", "<< xsec <<std::endl; 
      p_hist_tgt->SetBinContent(i_e,i_v,xsec*p_hist_tgt->GetBinContent(i_e,i_v)/total);

      //std::cout << p_hist->GetBinContent(i_e,i_v) << std::endl;
    }
  }
}

void CalcGENIECrossSections(){

  bool cc_or_nc = false;
  std::string rootdir = "/exp/uboone/data/users/jburridg/Nuance/";

  // Load the total cross section spline
  
  TFile* p_fxsec = TFile::Open("gxsec.root");
  TGraph* p_tot =  cc_or_nc ? static_cast<TGraph*>(p_fxsec->Get("nu_mu_C12/tot_cc")) : static_cast<TGraph*>(p_fxsec->Get("nu_mu_C12/tot_nc"));
  TGraph* p_tot_p =  cc_or_nc ? static_cast<TGraph*>(p_fxsec->Get("nu_mu_C12/tot_cc_p")) : static_cast<TGraph*>(p_fxsec->Get("nu_mu_C12/tot_nc_p"));
  TGraph* p_tot_n =  cc_or_nc ? static_cast<TGraph*>(p_fxsec->Get("nu_mu_C12/tot_cc_n")) : static_cast<TGraph*>(p_fxsec->Get("nu_mu_C12/tot_nc_n"));
  TGraph* p_tot_mec_pp = static_cast<TGraph*>(p_fxsec->Get("nu_mu_C12/mec_nc_pp")); //only applicable to NC interactions, no pp equivelnt for CC. 
  TGraph* p_tot_mec_nn =  cc_or_nc ? static_cast<TGraph*>(p_fxsec->Get("nu_mu_C12/mec_cc_nn")) : static_cast<TGraph*>(p_fxsec->Get("nu_mu_C12/mec_nc_nn"));
  TGraph* p_tot_mec_np =  cc_or_nc ? static_cast<TGraph*>(p_fxsec->Get("nu_mu_C12/mec_cc_np")) : static_cast<TGraph*>(p_fxsec->Get("nu_mu_C12/mec_nc_np"));
  
  
  // Convert from 10^-38 cm^2 to 10^-36 cm^2 as this is what NUANCE uses
  for(int i=0;i<p_tot->GetN();i++) p_tot->GetY()[i] /= 100; 
  gROOT->cd();
  p_fxsec->Close();

  std::string file = rootdir + "GENIE/" + (cc_or_nc ? "cc_events_low_e.root" : "nc_events_low_e.root");

  TFile* p_fin = TFile::Open(file.c_str());
  TTree* p_tin = static_cast<TTree*>(p_fin->Get("gst"));

  const int c_MAXPART = 25;

  // Setup branches
  Int_t           iev;
  Int_t           neu;
  Int_t           fspl;
  Int_t           tgt;
  Int_t           Z;
  Int_t           A;
  Int_t           hitnuc;
  Int_t           hitqrk;
  Int_t           resid;
  Bool_t          sea;
  Bool_t          qel;
  Bool_t          mec;
  Bool_t          res;
  Bool_t          dis;
  Bool_t          coh;
  Bool_t          dfr;
  Bool_t          imd;
  Bool_t          imdanh;
  Bool_t          singlek;
  Bool_t          nuel;
  Bool_t          em;
  Bool_t          cc;
  Bool_t          nc;
  Bool_t          charm;
  Bool_t          amnugamma;
  Bool_t          hnl;
  Int_t           neut_code;
  Int_t           nuance_code;
  Double_t        wght;
  Double_t        xs;
  Double_t        ys;
  Double_t        ts;
  Double_t        Q2s;
  Double_t        Ws;
  Double_t        x;
  Double_t        y;
  Double_t        t;
  Double_t        Q2;
  Double_t        W;
  Double_t        EvRF;
  Double_t        Ev;
  Double_t        pxv;
  Double_t        pyv;
  Double_t        pzv;
  Double_t        En;
  Double_t        pxn;
  Double_t        pyn;
  Double_t        pzn;
  Double_t        El;
  Double_t        pxl;
  Double_t        pyl;
  Double_t        pzl;
  Double_t        pl;
  Double_t        cthl;
  Int_t           nfp;
  Int_t           nfn;
  Int_t           nfpip;
  Int_t           nfpim;
  Int_t           nfpi0;
  Int_t           nfkp;
  Int_t           nfkm;
  Int_t           nfk0;
  Int_t           nfem;
  Int_t           nfother;
  Int_t           nip;
  Int_t           nin;
  Int_t           nipip;
  Int_t           nipim;
  Int_t           nipi0;
  Int_t           nikp;
  Int_t           nikm;
  Int_t           nik0;
  Int_t           niem;
  Int_t           niother;
  Int_t           ni;
  Int_t           pdgi[c_MAXPART];   //[ni]
  Int_t           resc[c_MAXPART];   //[ni]
  Double_t        Ei[c_MAXPART];   //[ni]
  Double_t        pxi[c_MAXPART];   //[ni]
  Double_t        pyi[c_MAXPART];   //[ni]
  Double_t        pzi[c_MAXPART];   //[ni]
  Int_t           nf;
  Int_t           pdgf[c_MAXPART];   //[nf]
  Double_t        Ef[c_MAXPART];   //[nf]
  Double_t        pxf[c_MAXPART];   //[nf]
  Double_t        pyf[c_MAXPART];   //[nf]
  Double_t        pzf[c_MAXPART];   //[nf]
  Double_t        pf[c_MAXPART];   //[nf]
  Double_t        cthf[c_MAXPART];   //[nf]
  Double_t        vtxx;
  Double_t        vtxy;
  Double_t        vtxz;
  Double_t        vtxt;
  Double_t        sumKEf;
  Double_t        calresp0;
  Double_t        XSec;
  Double_t        DXSec;
  UInt_t          KPS;
  p_tin->SetBranchAddress("iev", &iev);
  p_tin->SetBranchAddress("neu", &neu);
  p_tin->SetBranchAddress("fspl", &fspl);
  p_tin->SetBranchAddress("tgt", &tgt);
  p_tin->SetBranchAddress("Z", &Z);
  p_tin->SetBranchAddress("A", &A);
  p_tin->SetBranchAddress("hitnuc", &hitnuc);
  p_tin->SetBranchAddress("hitqrk", &hitqrk);
  p_tin->SetBranchAddress("resid", &resid);
  p_tin->SetBranchAddress("sea", &sea);
  p_tin->SetBranchAddress("qel", &qel);
  p_tin->SetBranchAddress("mec", &mec);
  p_tin->SetBranchAddress("res", &res);
  p_tin->SetBranchAddress("dis", &dis);
  p_tin->SetBranchAddress("coh", &coh);
  p_tin->SetBranchAddress("dfr", &dfr);
  p_tin->SetBranchAddress("imd", &imd);
  p_tin->SetBranchAddress("imdanh", &imdanh);
  p_tin->SetBranchAddress("singlek", &singlek);
  p_tin->SetBranchAddress("nuel", &nuel);
  p_tin->SetBranchAddress("em", &em);
  p_tin->SetBranchAddress("cc", &cc);
  p_tin->SetBranchAddress("nc", &nc);
  p_tin->SetBranchAddress("charm", &charm);
  p_tin->SetBranchAddress("amnugamma", &amnugamma);
  p_tin->SetBranchAddress("hnl", &hnl);
  p_tin->SetBranchAddress("neut_code", &neut_code);
  p_tin->SetBranchAddress("nuance_code", &nuance_code);
  p_tin->SetBranchAddress("wght", &wght);
  p_tin->SetBranchAddress("xs", &xs);
  p_tin->SetBranchAddress("ys", &ys);
  p_tin->SetBranchAddress("ts", &ts);
  p_tin->SetBranchAddress("Q2s", &Q2s);
  p_tin->SetBranchAddress("Ws", &Ws);
  p_tin->SetBranchAddress("x", &x);
  p_tin->SetBranchAddress("y", &y);
  p_tin->SetBranchAddress("t", &t);
  p_tin->SetBranchAddress("Q2", &Q2);
  p_tin->SetBranchAddress("W", &W);
  p_tin->SetBranchAddress("EvRF", &EvRF);
  p_tin->SetBranchAddress("Ev", &Ev);
  p_tin->SetBranchAddress("pxv", &pxv);
  p_tin->SetBranchAddress("pyv", &pyv);
  p_tin->SetBranchAddress("pzv", &pzv);
  p_tin->SetBranchAddress("En", &En);
  p_tin->SetBranchAddress("pxn", &pxn);
  p_tin->SetBranchAddress("pyn", &pyn);
  p_tin->SetBranchAddress("pzn", &pzn);
  p_tin->SetBranchAddress("El", &El);
  p_tin->SetBranchAddress("pxl", &pxl);
  p_tin->SetBranchAddress("pyl", &pyl);
  p_tin->SetBranchAddress("pzl", &pzl);
  p_tin->SetBranchAddress("pl", &pl);
  p_tin->SetBranchAddress("cthl", &cthl);
  p_tin->SetBranchAddress("nfp", &nfp);
  p_tin->SetBranchAddress("nfn", &nfn);
  p_tin->SetBranchAddress("nfpip", &nfpip);
  p_tin->SetBranchAddress("nfpim", &nfpim);
  p_tin->SetBranchAddress("nfpi0", &nfpi0);
  p_tin->SetBranchAddress("nfkp", &nfkp);
  p_tin->SetBranchAddress("nfkm", &nfkm);
  p_tin->SetBranchAddress("nfk0", &nfk0);
  p_tin->SetBranchAddress("nfem", &nfem);
  p_tin->SetBranchAddress("nfother", &nfother);
  p_tin->SetBranchAddress("nip", &nip);
  p_tin->SetBranchAddress("nin", &nin);
  p_tin->SetBranchAddress("nipip", &nipip);
  p_tin->SetBranchAddress("nipim", &nipim);
  p_tin->SetBranchAddress("nipi0", &nipi0);
  p_tin->SetBranchAddress("nikp", &nikp);
  p_tin->SetBranchAddress("nikm", &nikm);
  p_tin->SetBranchAddress("nik0", &nik0);
  p_tin->SetBranchAddress("niem", &niem);
  p_tin->SetBranchAddress("niother", &niother);
  p_tin->SetBranchAddress("ni", &ni);
  p_tin->SetBranchAddress("pdgi", pdgi);
  p_tin->SetBranchAddress("resc", resc);
  p_tin->SetBranchAddress("Ei", Ei);
  p_tin->SetBranchAddress("pxi", pxi);
  p_tin->SetBranchAddress("pyi", pyi);
  p_tin->SetBranchAddress("pzi", pzi);
  p_tin->SetBranchAddress("nf", &nf);
  p_tin->SetBranchAddress("pdgf", pdgf);
  p_tin->SetBranchAddress("Ef", Ef);
  p_tin->SetBranchAddress("pxf", pxf);
  p_tin->SetBranchAddress("pyf", pyf);
  p_tin->SetBranchAddress("pzf", pzf);
  p_tin->SetBranchAddress("pf", pf);
  p_tin->SetBranchAddress("cthf", cthf);
  p_tin->SetBranchAddress("vtxx", &vtxx);
  p_tin->SetBranchAddress("vtxy", &vtxy);
  p_tin->SetBranchAddress("vtxz", &vtxz);
  p_tin->SetBranchAddress("vtxt", &vtxt);
  p_tin->SetBranchAddress("sumKEf", &sumKEf);
  p_tin->SetBranchAddress("calresp0", &calresp0);
  p_tin->SetBranchAddress("XSec", &XSec);
  p_tin->SetBranchAddress("DXSec", &DXSec);
  p_tin->SetBranchAddress("KPS", &KPS);

  const Long64_t c_nevents = p_tin->GetEntries();
  //std::cout << "Tree has " << c_nevents << " events" << std::endl;
  const int spline_range = 2.0 //helps keep all spline referencing consistent
  // Record the total number of events generated for each neutrino energy
  TH1D* p_nevents = new TH1D("nevents",";Neutrino Energy (GeV);N Events",40,0.0,10.0);
  TH1D* p_nevents_p = new TH1D("nevents_p",";Neutrino Energy (Proton target) (GeV);N Events",40,0.0,spline_range);
  TH1D* p_nevents_n = new TH1D("nevents_n",";Neutrino Energy (Neutron target) (GeV);N Events",40,0.0,spline_range);

  TH1D* p_nevents_mec_nn = new TH1D("nevents_mec_nn",";Neutrino Energy (MEC nn)(GeV);N Events",40,0.0,spline_range);
  TH1D* p_nevents_mec_pp = new TH1D("nevents_mec_pp",";Neutrino Energy (MEC pp) (GeV);N Events",40,0.0,spline_range);
  TH1D* p_nevents_mec_np = new TH1D("nevents_mec_np",";Neutrino Energy (MEC np) (GeV);N Events",40,0.0,spline_range);


  // Setup histograms - to normalise using the xsec spline we always set the x axis to be 
  // in bins of neutrino energy and the y axis as whatever other varuable we want
  TH2D* p_leptonmomentum_p = new TH2D("muonmomentum","Inclusive (Proton);Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
  TH2D* p_leptonmomentum_n = new TH2D("muonmomentum","Inclusive (Neutron);Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
  TH2D* p_leptoncostheta_p = new TH2D("muoncostheta","Inclusive (Proton);Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);
  TH2D* p_leptoncostheta_n = new TH2D("muoncostheta","Inclusive (Neutron);Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);

  //MEC Interactions:
  TH2D* p_leptonmomentum_mec_nn = new TH2D("muonmomentum","MEC (nn);Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
  TH2D* p_leptonmomentum_mec_np = new TH2D("muonmomentum","MEC (np);Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
  TH2D* p_leptonmomentum_mec_pp = new TH2D("muonmomentum","MEC (pp);Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
  TH2D* p_leptoncostheta_mec_nn = new TH2D("muoncostheta","MEC (nn);Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);
  TH2D* p_leptoncostheta_mec_np = new TH2D("muoncostheta","MEC (np);Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);
  TH2D* p_leptoncostheta_mec_pp = new TH2D("muoncostheta","MEC (pp);Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);
  

  std::map<int,TH2D*> m_ch_leptonmomentum;
  std::map<int,TH2D*> m_ch_leptonmomentum_p;
  std::map<int,TH2D*> m_ch_leptonmomentum_n;
  std::map<int,TH2D*> m_ch_leptoncostheta_p;
  std::map<int,TH2D*> m_ch_leptoncostheta_n;

  for(Long64_t ievent=0;ievent<c_nevents;ievent++){
    p_tin->GetEntry(ievent);
    //std::cout<< Ev << std::endl;
    //std::cout << nuance_code << std::endl;
    //std::cout<< mec << "  " << hitnuc << " " << std::endl; 

    p_nevents->Fill(Ev);

    double mom = sqrt(pxl*pxl+pyl*pyl+pzl*pzl);
    double costheta = pzl/mom;
    //std::cout<< pdgi <<std::endl;

    if (hitnuc == 2212){
      p_nevents_p->Fill(Ev);
      p_leptonmomentum_p->Fill(Ev,mom);
      p_leptoncostheta_p->Fill(Ev,costheta);

      if(m_ch_leptonmomentum_p.find(nuance_code) == m_ch_leptonmomentum_p.end()){
        m_ch_leptonmomentum_p[nuance_code] = new TH2D(Form("muonmomentum_p_%i",nuance_code),";Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
        m_ch_leptoncostheta_p[nuance_code] = new TH2D(Form("muoncostheta_p_%i",nuance_code),";Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);  
      }
      m_ch_leptonmomentum_p.at(nuance_code)->Fill(Ev,mom);
      m_ch_leptoncostheta_p.at(nuance_code)->Fill(Ev,costheta);

    } else if (hitnuc == 2112){
      
      p_nevents_n->Fill(Ev);
      p_leptonmomentum_n->Fill(Ev,mom);
      p_leptoncostheta_n->Fill(Ev,costheta);

      if(m_ch_leptonmomentum_n.find(nuance_code) == m_ch_leptonmomentum_n.end()){
        m_ch_leptonmomentum_n[nuance_code] = new TH2D(Form("muonmomentum_n_%i",nuance_code),";Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
        m_ch_leptoncostheta_n[nuance_code] = new TH2D(Form("muoncostheta_n_%i",nuance_code),";Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);
      }
      
      m_ch_leptonmomentum_n.at(nuance_code)->Fill(Ev,mom);
      m_ch_leptoncostheta_n.at(nuance_code)->Fill(Ev,costheta);
      
    } 
    
    if (cc_or_nc){
      // CC
      if (hitnuc == 2000000200){ //nn
        p_nevents_mec_nn->Fill(Ev);
        p_leptonmomentum_mec_nn->Fill(Ev,mom);
        p_leptoncostheta_mec_nn->Fill(Ev,costheta);
      }
      else if (hitnuc == 2000000201){ //np
        p_nevents_mec_np->Fill(Ev);
        p_leptonmomentum_mec_np->Fill(Ev,mom);
        p_leptoncostheta_mec_np->Fill(Ev,costheta);
      }
      
     // std::cout<< mec << "  " << hitnuc << " " << std::endl; 
    } else{
      //NC
      if (hitnuc == 2000000200){ //nn
        p_nevents_mec_nn->Fill(Ev);
        p_leptonmomentum_mec_nn->Fill(Ev,mom);
        p_leptoncostheta_mec_nn->Fill(Ev,costheta);
      }
      else if (hitnuc == 2000000201){ //np
        p_nevents_mec_np->Fill(Ev);
        p_leptonmomentum_mec_np->Fill(Ev,mom);
        p_leptoncostheta_mec_np->Fill(Ev,costheta);
      }
      else if (hitnuc == 2000000202){ //pp
        p_nevents_mec_pp->Fill(Ev);
        p_leptonmomentum_mec_pp->Fill(Ev,mom);
        p_leptoncostheta_mec_pp->Fill(Ev,costheta);
      }
    }
  }

  //NormaliseToSpline(p_tot_p,p_nevents_p,p_leptonmomentum_p);
  //NormaliseToSpline(p_tot_n,p_nevents_n,p_leptonmomentum_n);
 // NormaliseToSpline(p_tot_p,p_nevents_p,p_leptoncostheta_p);
  //NormaliseToSpline(p_tot_n,p_nevents_n,p_leptoncostheta_n);

  //maybe add a conditional here about if its NC or CC? need to know weather to fill mec_pp or not. 
  NormaliseToSpline(p_tot_mec_pp,p_nevents_mec_pp,p_leptonmomentum_mec_pp);
  NormaliseToSpline(p_tot_mec_np,p_nevents_mec_np,p_leptonmomentum_mec_np);
  NormaliseToSpline(p_tot_mec_nn,p_nevents_mec_nn,p_leptonmomentum_mec_nn);
  NormaliseToSpline(p_tot_mec_pp,p_nevents_mec_pp,p_leptoncostheta_mec_pp);
  NormaliseToSpline(p_tot_mec_np,p_nevents_mec_np,p_leptoncostheta_mec_np);
  NormaliseToSpline(p_tot_mec_nn,p_nevents_mec_nn,p_leptoncostheta_mec_nn);
 


  TFile* p_fout = cc_or_nc ? new TFile("GENIECrossSections_CC.root","RECREATE") : new TFile("GENIECrossSections_NC.root","RECREATE");
  p_fout->cd();

  for(std::map<int,TH2D*>::iterator it = m_ch_leptonmomentum_p.begin();it != m_ch_leptonmomentum_p.end();it++){
    NormaliseToSpline(p_tot_p,p_nevents_p,it->second);
    it->second->Write(Form("leptonmomentum_p_%i",it->first));
  }
  for(std::map<int,TH2D*>::iterator it = m_ch_leptonmomentum_n.begin();it != m_ch_leptonmomentum_n.end();it++){
    NormaliseToSpline(p_tot_n,p_nevents_n,it->second);
    it->second->Write(Form("leptonmomentum_n_%i",it->first));
  }
  for(std::map<int,TH2D*>::iterator it = m_ch_leptoncostheta_p.begin();it != m_ch_leptoncostheta_p.end();it++){
    NormaliseToSpline(p_tot_p,p_nevents_p,it->second);
    it->second->Write(Form("leptoncostheta_p_%i",it->first));
  }
  for(std::map<int,TH2D*>::iterator it = m_ch_leptoncostheta_n.begin();it != m_ch_leptoncostheta_n.end();it++){
    NormaliseToSpline(p_tot_n,p_nevents_n,it->second);
    it->second->Write(Form("leptoncostheta_n_%i",it->first));
  }

  if(cc_or_nc) {
    p_leptonmomentum_p->Write("leptonmomentum_cc_p");
    p_leptonmomentum_n->Write("leptonmomentum_cc_n");
    p_leptonmomentum_mec_np->Write("leptonmomentum_cc_mec_np");
    p_leptonmomentum_mec_nn->Write("leptonmomentum_cc_mec_nn");
  }
  else {
    p_leptonmomentum_p->Write("leptonmomentum_nc_p");
    p_leptonmomentum_n->Write("leptonmomentum_nc_n");
    p_leptonmomentum_mec_np->Write("leptonmomentum_nc_mec_np");
    p_leptonmomentum_mec_nn->Write("leptonmomentum_nc_mec_nn");
    p_leptonmomentum_mec_pp->Write("leptonmomentum_nc_mec_pp");
  }
  if(cc_or_nc) {
    p_leptoncostheta_p->Write("leptoncostheta_cc_p");
    p_leptoncostheta_n->Write("leptoncostheta_cc_n");
    p_leptoncostheta_mec_np->Write("leptoncostheta_cc_mec_np");
    p_leptoncostheta_mec_nn->Write("leptoncostheta_cc_mec_nn");
  }
  else {
    p_leptoncostheta_p->Write("leptoncostheta_nc_p");
    p_leptoncostheta_n->Write("leptoncostheta_nc_n");
    p_leptoncostheta_mec_np->Write("leptoncostheta_nc_mec_np");
    p_leptoncostheta_mec_nn->Write("leptoncostheta_nc_mec_nn");
    p_leptoncostheta_mec_pp->Write("leptoncostheta_nc_mec_pp");
  }

  TCanvas* p_canvas = new TCanvas("c","c");

  p_leptonmomentum_p->Draw("colz");
  p_leptonmomentum_p->SetStats(0);
  if(cc_or_nc) p_canvas->Print("GENIE_C_leptonmomentum_cc_p.png");
  else p_canvas->Print("GENIE_C_leptonmomentum_nc_p.png");
  p_canvas->Clear();  

  p_leptonmomentum_n->Draw("colz");
  p_leptonmomentum_n->SetStats(0);
  if(cc_or_nc) p_canvas->Print("GENIE_C_leptonmomentum_cc_n.png");
  else p_canvas->Print("GENIE_C_leptonmomentum_nc_n.png");
  p_canvas->Clear(); 

  p_leptoncostheta_p->Draw("colz");
  p_leptoncostheta_p->SetStats(0);
  if(cc_or_nc) p_canvas->Print("GENIE_C_leptoncostheta_cc_p.png");
  else p_canvas->Print("GENIE_C_leptoncostheta_nc_p.png");
  p_canvas->Clear();  

  p_leptoncostheta_n->Draw("colz");
  p_leptoncostheta_n->SetStats(0);
  if(cc_or_nc) p_canvas->Print("GENIE_C_leptoncostheta_cc_n.png");
  else p_canvas->Print("GENIE_C_leptoncostheta_nc_n.png");
  p_canvas->Clear();  

  //MEC Histograms
  p_leptonmomentum_mec_pp->Draw("colz"); //no CC pp interaction. 
  p_leptonmomentum_mec_pp->SetStats(0);
  if(!cc_or_nc) p_canvas->Print("GENIE_C_leptonmomentum_nc_mec_pp.png");
  p_canvas->Clear();  

  p_leptoncostheta_mec_pp->Draw("colz");
  p_leptoncostheta_mec_pp->SetStats(0);
  if(!cc_or_nc) p_canvas->Print("GENIE_C_leptoncostheta_nc_mec_pp.png");
  p_canvas->Clear();  

  p_leptonmomentum_mec_np->Draw("colz");
  p_leptonmomentum_mec_np->SetStats(0);
  if(cc_or_nc) p_canvas->Print("GENIE_C_leptonmomentum_cc_mec_np.png");
  else p_canvas->Print("GENIE_C_leptonmomentum_nc_mec_np.png");
  p_canvas->Clear(); 
  
  p_leptoncostheta_mec_np->Draw("colz");
  p_leptoncostheta_mec_np->SetStats(0);
  if(cc_or_nc) p_canvas->Print("GENIE_C_leptoncostheta_cc_mec_np.png");
  else p_canvas->Print("GENIE_C_leptoncostheta_nc_mec_np.png");
  p_canvas->Clear();  

  p_leptonmomentum_mec_nn->Draw("colz");
  p_leptonmomentum_mec_nn->SetStats(0);
  if(cc_or_nc) p_canvas->Print("GENIE_C_leptonmomentum_cc_mec_nn.png");
  else p_canvas->Print("GENIE_C_leptonmomentum_nc_mec_nn.png");
  p_canvas->Clear();  

  p_leptoncostheta_mec_nn->Draw("colz");
  p_leptoncostheta_mec_nn->SetStats(0);
  if(cc_or_nc) p_canvas->Print("GENIE_C_leptoncostheta_cc_mec_nn.png");
  else p_canvas->Print("GENIE_C_leptoncostheta_nc_mec_nn.png");
  p_canvas->Clear();  
  
  //1D Histograms
  p_nevents_p->Draw();
  p_canvas->Print("GENIE_C_pnevents_p.png");
  p_canvas->Clear(); 

  p_nevents_n->Draw();
  p_canvas->Print("GENIE_C_pnevents_n.png");
  p_canvas->Clear();

  p_nevents->Draw();
  p_canvas->Print("GENIE_C_pnevents.png");
  p_canvas->Clear(); 
  
  //MEC 1D Histograms
  p_nevents_mec_pp->Draw();
  if(!cc_or_nc) p_canvas->Print("GENIE_C_pnevents_nc_mec_pp.png");
  p_canvas->Clear(); 

  p_nevents_mec_nn->Draw();
  if(cc_or_nc) p_canvas->Print("GENIE_C_pnevents_cc_mec_nn.png");
  else p_canvas->Print("GENIE_C_pnevents_nc_mec_nn.png");
  p_canvas->Clear();  

  p_nevents_mec_np->Draw();
  if(cc_or_nc) p_canvas->Print("GENIE_C_pnevents_cc_mec_np.png");
  else p_canvas->Print("GENIE_C_pnevents_nc_mec_np.png");
  p_canvas->Clear();  


  p_fout->Close();
  p_fin->Close();
  
  
}
