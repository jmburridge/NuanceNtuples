
void DrawTotalCrossSections(){

  // First extract the h100 tree from the file

  TFile* p_fin = TFile::Open("nuance_v3_may07_20070507_numu_ma1.35_kappa1.007_1.root");
  TTree* p_tin = static_cast<TTree*>(p_fin->Get("h100"));

  Float_t         Enu;
  Float_t         num_pfcc;
  Float_t         num_pfnc;
  Float_t         num_ecc;
  Float_t         num_enc;
  Float_t         num_ccc;
  Float_t         num_cnc;
  Float_t         num_nbcc;
  Float_t         num_nbnc;
  Float_t         num_pbcc;
  Float_t         num_pbnc;
  p_tin->SetBranchAddress("Enu", &Enu);
  p_tin->SetBranchAddress("num_pfcc", &num_pfcc);
  p_tin->SetBranchAddress("num_pfnc", &num_pfnc);
  p_tin->SetBranchAddress("num_ecc", &num_ecc);
  p_tin->SetBranchAddress("num_enc", &num_enc);
  p_tin->SetBranchAddress("num_ccc", &num_ccc);
  p_tin->SetBranchAddress("num_cnc", &num_cnc);
  p_tin->SetBranchAddress("num_nbcc", &num_nbcc);
  p_tin->SetBranchAddress("num_nbnc", &num_nbnc);
  p_tin->SetBranchAddress("num_pbcc", &num_pbcc);
  p_tin->SetBranchAddress("num_pbnc", &num_pbnc);

  std::vector<Float_t> v_Enu; 
  std::vector<Float_t> v_num_pfcc; 
  std::vector<Float_t> v_num_pfnc; 
  std::vector<Float_t> v_num_ecc; 
  std::vector<Float_t> v_num_enc; 
  std::vector<Float_t> v_num_ccc; 
  std::vector<Float_t> v_num_cnc; 
  std::vector<Float_t> v_num_nbcc; 
  std::vector<Float_t> v_num_nbnc; 
  std::vector<Float_t> v_num_pbcc; 
  std::vector<Float_t> v_num_pbnc; 

  // Couple of extra vectors to store the C xsec/12 for the plots and the total C xsec
  std::vector<Float_t> v_num_c; 
  std::vector<Float_t> v_num_ccc_12; 
  std::vector<Float_t> v_num_cnc_12; 

  // Convert the data from the tree into TGraphs and save
  const Long64_t c_nentries = p_tin->GetEntries();

  for(Long64_t ientry=0;ientry<c_nentries;ientry++){
    p_tin->GetEntry(ientry);
    v_Enu.push_back(Enu);      
    v_num_pfcc.push_back(num_pfcc);      
    v_num_pfnc.push_back(num_pfnc);      
    v_num_ecc.push_back(num_ecc);      
    v_num_enc.push_back(num_enc);      
    v_num_ccc.push_back(num_ccc);      
    v_num_cnc.push_back(num_cnc);
    v_num_nbcc.push_back(num_nbcc);
    v_num_nbnc.push_back(num_nbnc);      
    v_num_pbcc.push_back(num_pbcc);      
    v_num_pbnc.push_back(num_pbnc);      
    v_num_c.push_back(num_ccc+num_cnc);      
    v_num_ccc_12.push_back(num_ccc/12);      
    v_num_cnc_12.push_back(num_cnc/12);
  }

  p_fin->Close(); 

  TFile* p_fout = new TFile("NuanceSplines.root","RECREATE");
  TGraph* p_g_num_pfcc = new TGraph(v_Enu.size(),&(v_Enu[0]),&(v_num_pfcc[0]));
  TGraph* p_g_num_pfnc = new TGraph(v_Enu.size(),&(v_Enu[0]),&(v_num_pfnc[0]));
  TGraph* p_g_num_ecc = new TGraph(v_Enu.size(),&(v_Enu[0]),&(v_num_ecc[0]));
  TGraph* p_g_num_enc = new TGraph(v_Enu.size(),&(v_Enu[0]),&(v_num_enc[0]));
  TGraph* p_g_num_ccc = new TGraph(v_Enu.size(),&(v_Enu[0]),&(v_num_ccc[0]));
  TGraph* p_g_num_cnc = new TGraph(v_Enu.size(),&(v_Enu[0]),&(v_num_cnc[0]));
  TGraph* p_g_num_nbcc = new TGraph(v_Enu.size(),&(v_Enu[0]),&(v_num_nbcc[0]));
  TGraph* p_g_num_nbnc = new TGraph(v_Enu.size(),&(v_Enu[0]),&(v_num_nbnc[0]));
  TGraph* p_g_num_pbcc = new TGraph(v_Enu.size(),&(v_Enu[0]),&(v_num_pbcc[0]));
  TGraph* p_g_num_pbnc = new TGraph(v_Enu.size(),&(v_Enu[0]),&(v_num_pbnc[0]));
  TGraph* p_g_num_c = new TGraph(v_Enu.size(),&(v_Enu[0]),&(v_num_c[0]));
  TGraph* p_g_num_ccc_12 = new TGraph(v_Enu.size(),&(v_Enu[0]),&(v_num_ccc_12[0]));
  TGraph* p_g_num_cnc_12 = new TGraph(v_Enu.size(),&(v_Enu[0]),&(v_num_cnc_12[0]));
   
  TMultiGraph* p_mg = new TMultiGraph("mg",";E_{#nu} (MeV);#sigma (10^{-36} cm^{2})");
  TLegend* p_leg = new TLegend(0.1,0.7,0.5,0.9); 
  p_leg->SetNColumns(2);
 
  p_g_num_pfcc->SetLineColor(1);
  p_g_num_pfcc->SetLineStyle(1);
  p_mg->Add(p_g_num_pfcc);
  p_leg->AddEntry(p_g_num_pfcc,"Per free nucleon","L");
  
  p_g_num_pfnc->SetLineColor(1);
  p_g_num_pfnc->SetLineStyle(2);
  p_mg->Add(p_g_num_pfnc);
   
  p_g_num_ecc->SetLineColor(2);
  p_g_num_ecc->SetLineStyle(1);
  p_mg->Add(p_g_num_ecc);
  p_leg->AddEntry(p_g_num_ecc,"Per electron","L");

  p_g_num_enc->SetLineColor(2);
  p_g_num_enc->SetLineStyle(2);
  p_mg->Add(p_g_num_enc);

  p_g_num_ccc_12->SetLineColor(3);
  p_g_num_ccc_12->SetLineStyle(1);
  p_mg->Add(p_g_num_ccc_12);
  p_leg->AddEntry(p_g_num_ccc_12,"Per C nucleus/12","L");

  p_g_num_cnc_12->SetLineColor(3);
  p_g_num_cnc_12->SetLineStyle(2);
  p_mg->Add(p_g_num_cnc_12);

  p_g_num_nbcc->SetLineColor(4);
  p_g_num_nbcc->SetLineStyle(1);
  p_mg->Add(p_g_num_nbcc);
  p_leg->AddEntry(p_g_num_nbcc,"Per bound neutron","L");

  p_g_num_nbnc->SetLineColor(4);
  p_g_num_nbnc->SetLineStyle(2);
  p_mg->Add(p_g_num_nbnc);

  p_g_num_pbcc->SetLineColor(5);
  p_g_num_pbcc->SetLineStyle(1);
  p_mg->Add(p_g_num_pbcc);
  p_leg->AddEntry(p_g_num_pbcc,"Per bound proton","L");

  p_g_num_pbnc->SetLineColor(5);
  p_g_num_pbnc->SetLineStyle(2);
  p_mg->Add(p_g_num_pbnc);

  TCanvas* p_can = new TCanvas("can","can");
  p_mg->Draw("AL"); 
  p_leg->Draw();
  p_can->Print("test.png");

  // Write everything to a spline file 
  p_g_num_pfcc->Write("num_pfcc");
  p_g_num_pfnc->Write("num_pfnc");
  p_g_num_ecc->Write("num_ecc");
  p_g_num_enc->Write("num_enc");
  p_g_num_ccc->Write("num_ccc");
  p_g_num_cnc->Write("num_cnc");
  p_g_num_nbcc->Write("num_nbcc");
  p_g_num_nbnc->Write("num_nbnc");
  p_g_num_pbcc->Write("num_pbcc");
  p_g_num_pbnc->Write("num_pbnc");
  p_g_num_c->Write("num_c");

  p_fout->Close();

}
