
void CalculateRatios(){

  // Load the NUANCE Cross Sections 
  TFile* p_fin_NUANCE = TFile::Open("NuanceCrossSections.root");
  TH2D* p_NUANCE_leptonmomentum_cc_p = static_cast<TH2D*>(p_fin_NUANCE->Get("leptonmomentum_cc_p"));
  TH2D* p_NUANCE_leptonmomentum_nc_p = static_cast<TH2D*>(p_fin_NUANCE->Get("leptonmomentum_nc_p"));
  TH2D* p_NUANCE_leptonmomentum_cc_n = static_cast<TH2D*>(p_fin_NUANCE->Get("leptonmomentum_cc_n"));
  TH2D* p_NUANCE_leptonmomentum_nc_n = static_cast<TH2D*>(p_fin_NUANCE->Get("leptonmomentum_nc_n"));

  TH2D* p_NUANCE_leptoncostheta_cc_p = static_cast<TH2D*>(p_fin_NUANCE->Get("leptoncostheta_cc_p"));
  TH2D* p_NUANCE_leptoncostheta_nc_p = static_cast<TH2D*>(p_fin_NUANCE->Get("leptoncostheta_nc_p"));
  TH2D* p_NUANCE_leptoncostheta_cc_n = static_cast<TH2D*>(p_fin_NUANCE->Get("leptoncostheta_cc_n"));
  TH2D* p_NUANCE_leptoncostheta_nc_n = static_cast<TH2D*>(p_fin_NUANCE->Get("leptoncostheta_nc_n"));

  p_NUANCE_leptonmomentum_cc_p->SetDirectory(0);
  p_NUANCE_leptonmomentum_nc_p->SetDirectory(0);
  p_NUANCE_leptonmomentum_cc_n->SetDirectory(0);
  p_NUANCE_leptonmomentum_nc_n->SetDirectory(0);

  p_NUANCE_leptoncostheta_cc_p->SetDirectory(0);
  p_NUANCE_leptoncostheta_nc_p->SetDirectory(0);
  p_NUANCE_leptoncostheta_cc_n->SetDirectory(0);
  p_NUANCE_leptoncostheta_nc_n->SetDirectory(0);
  p_fin_NUANCE->Close();


  // Load GENIE Cross Sections
  TFile* p_fin_GENIE_cc = TFile::Open("GENIECrossSections_CC.root");
  TH2D* p_GENIE_leptonmomentum_cc_p = static_cast<TH2D*>(p_fin_GENIE_cc->Get("leptonmomentum_cc_p"));
  TH2D* p_GENIE_leptonmomentum_cc_n = static_cast<TH2D*>(p_fin_GENIE_cc->Get("leptonmomentum_cc_n"));
  TH2D* p_GENIE_leptoncostheta_cc_p = static_cast<TH2D*>(p_fin_GENIE_cc->Get("leptoncostheta_cc_p"));
  TH2D* p_GENIE_leptoncostheta_cc_n = static_cast<TH2D*>(p_fin_GENIE_cc->Get("leptoncostheta_cc_n"));
  p_GENIE_leptonmomentum_cc_p->SetDirectory(0);
  p_GENIE_leptonmomentum_cc_n->SetDirectory(0);
  p_GENIE_leptoncostheta_cc_p->SetDirectory(0);
  p_GENIE_leptoncostheta_cc_n->SetDirectory(0);
  p_fin_GENIE_cc->Close();

  TFile* p_fin_GENIE_nc = TFile::Open("GENIECrossSections_NC.root");
  TH2D* p_GENIE_leptonmomentum_nc_p = static_cast<TH2D*>(p_fin_GENIE_nc->Get("leptonmomentum_nc_p"));
  TH2D* p_GENIE_leptonmomentum_nc_n = static_cast<TH2D*>(p_fin_GENIE_nc->Get("leptonmomentum_nc_n"));
  TH2D* p_GENIE_leptoncostheta_nc_p = static_cast<TH2D*>(p_fin_GENIE_nc->Get("leptoncostheta_nc_p"));
  TH2D* p_GENIE_leptoncostheta_nc_n= static_cast<TH2D*>(p_fin_GENIE_nc->Get("leptoncostheta_nc_n"));
  p_GENIE_leptonmomentum_nc_p->SetDirectory(0);
  p_GENIE_leptonmomentum_nc_n->SetDirectory(0);
  p_GENIE_leptoncostheta_nc_p->SetDirectory(0);
  p_GENIE_leptoncostheta_nc_n->SetDirectory(0);
  p_fin_GENIE_nc->Close();

  p_NUANCE_leptonmomentum_cc_p->Divide(p_GENIE_leptonmomentum_cc_p);
  p_NUANCE_leptonmomentum_nc_p->Divide(p_GENIE_leptonmomentum_nc_p);
  p_NUANCE_leptonmomentum_cc_n->Divide(p_GENIE_leptonmomentum_cc_n);
  p_NUANCE_leptonmomentum_nc_n->Divide(p_GENIE_leptonmomentum_nc_n);

  p_NUANCE_leptoncostheta_cc_p->Divide(p_GENIE_leptoncostheta_cc_p);
  p_NUANCE_leptoncostheta_nc_p->Divide(p_GENIE_leptoncostheta_nc_p);
  p_NUANCE_leptoncostheta_cc_n->Divide(p_GENIE_leptoncostheta_cc_n);
  p_NUANCE_leptoncostheta_nc_n->Divide(p_GENIE_leptoncostheta_nc_n);
  
  TFile* p_ratios = new TFile("Ratios.root","RECREATE");
  p_ratios->cd();  
  p_NUANCE_leptonmomentum_cc_p->Write("leptonmomentum_cc_p"); 
  p_NUANCE_leptonmomentum_nc_p->Write("leptonmomentum_nc_p"); 
  p_NUANCE_leptonmomentum_cc_n->Write("leptonmomentum_cc_n"); 
  p_NUANCE_leptonmomentum_nc_n->Write("leptonmomentum_nc_n"); 

  p_NUANCE_leptoncostheta_cc_p->Write("leptoncostheta_cc_p"); 
  p_NUANCE_leptoncostheta_nc_p->Write("leptoncostheta_nc_p"); 
  p_NUANCE_leptoncostheta_cc_n->Write("leptoncostheta_cc_n"); 
  p_NUANCE_leptoncostheta_nc_n->Write("leptoncostheta_nc_n"); 
  p_ratios->Close();

}
