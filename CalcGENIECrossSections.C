#include "TFile.h"
#include "TTree.h"

// Normalise each vertical strip of a 2D histogram 
// Inputs: p_spline gives the total cross section for C nuclei, p_total is the total number of events
// on C nuclei in the ntuple for each neutrino energy, p_hist is the distribution of events in the ntuple in
// neutrino energy and some other variable on the y axis
void NormaliseToSpline(TGraph* p_spline,TH1D* p_total,TH2D* p_hist){
  //std::cout << "Normalising to spline" << std::endl;
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

void CalcGENIECrossSections(){

  bool cc_or_nc = false;

  // Load the total cross section spline
  
  TFile* p_fxsec = TFile::Open("gxsec.root");
  TGraph* p_tot =  cc_or_nc ? static_cast<TGraph*>(p_fxsec->Get("nu_mu_C12/tot_cc")) : static_cast<TGraph*>(p_fxsec->Get("nu_mu_C12/tot_nc"));
  // Convert from 10^-38 cm^2 to 10^-36 cm^2 as this is what NUANCE uses
  for(int i=0;i<p_tot->GetN();i++) p_tot->GetY()[i] /= 100; 
  gROOT->cd();
  p_fxsec->Close();

  TFile* p_fin = cc_or_nc ? TFile::Open("cc_events.root") : TFile::Open("nc_events.root");
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

  // Record the total number of events generated for each neutrino energy
  TH1D* p_nevents = new TH1D("nevents",";Neutrino Energy (GeV);N Events",40,0.0,2.0);

  // Setup histograms - to normalise using the xsec spline we always set the x axis to be 
  // in bins of neutrino energy and the y axis as whatever other varuable we want
  TH2D* p_leptonmomentum = new TH2D("muonmomentum","Inclusive;Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
  TH2D* p_leptoncostheta = new TH2D("muoncostheta","Inclusive;Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);

  std::map<int,TH2D*> m_ch_leptonmomentum;
  std::map<int,TH2D*> m_ch_leptoncostheta;

  for(Long64_t ievent=0;ievent<c_nevents;ievent++){
    p_tin->GetEntry(ievent);
    
    //std::cout << nuance_code << std::endl;

    p_nevents->Fill(Ev);

    double mom = sqrt(pxl*pxl+pyl*pyl+pzl*pzl);
    p_leptonmomentum->Fill(Ev,mom);

    double costheta = pzl/mom;
    p_leptoncostheta->Fill(Ev,costheta);

    if(m_ch_leptonmomentum.find(nuance_code) == m_ch_leptonmomentum.end()){
      m_ch_leptonmomentum[nuance_code] = new TH2D(Form("muonmomentum_%i",nuance_code),";Neutrino Energy (GeV);Lepton Momentum (GeV);d#sigma/dP (10^{-36} cm^2/GeV)",40,0.0,2.0,40,0.0,2.0);
      m_ch_leptoncostheta[nuance_code] = new TH2D(Form("muoncostheta_%i",nuance_code),";Neutrino Energy (GeV);Lepton Cos(#theta);d#sigma/dCos(#theta) (10^{-36} cm^2)",40,0.0,2.0,40,-1.0,1.0);
    }

    m_ch_leptonmomentum.at(nuance_code)->Fill(Ev,mom);
    m_ch_leptoncostheta.at(nuance_code)->Fill(Ev,costheta);
      
  }

  NormaliseToSpline(p_tot,p_nevents,p_leptonmomentum);
  NormaliseToSpline(p_tot,p_nevents,p_leptoncostheta);

  TFile* p_fout = cc_or_nc ? new TFile("GENIECrossSections_CC.root","RECREATE") : new TFile("GENIECrossSections_NC.root","RECREATE");
  p_fout->cd();

  for(std::map<int,TH2D*>::iterator it = m_ch_leptonmomentum.begin();it != m_ch_leptonmomentum.end();it++){
    NormaliseToSpline(p_tot,p_nevents,it->second);
    it->second->Write(Form("leptonmomentum_%i",it->first));
  }
  for(std::map<int,TH2D*>::iterator it = m_ch_leptoncostheta.begin();it != m_ch_leptoncostheta.end();it++){
    NormaliseToSpline(p_tot,p_nevents,it->second);
    it->second->Write(Form("leptoncostheta_%i",it->first));
  }

  if(cc_or_nc) p_leptonmomentum->Write("leptonmomentum_cc");
  else p_leptonmomentum->Write("leptonmomentum_nc");
  if(cc_or_nc) p_leptoncostheta->Write("leptoncostheta_cc");
  else p_leptoncostheta->Write("leptoncostheta_nc");

  TCanvas* p_canvas = new TCanvas("c","c");

  p_leptonmomentum->Draw("colz");
  p_leptonmomentum->SetStats(0);
  if(cc_or_nc) p_canvas->Print("GENIE_leptonmomentum_cc.png");
  else p_canvas->Print("GENIE_leptonmomentum_nc.png");
  p_canvas->Clear();  

  p_leptoncostheta->Draw("colz");
  p_leptoncostheta->SetStats(0);
  if(cc_or_nc) p_canvas->Print("GENIE_leptoncostheta_cc.png");
  else p_canvas->Print("GENIE_leptoncostheta_nc.png");
  p_canvas->Clear();  

  p_fout->Close();
  p_fin->Close();
  
}
