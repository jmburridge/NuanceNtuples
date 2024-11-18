
void CalculateRatios(){

  // Load the NUANCE Cross Sections 
  TFile* p_fin_NUANCE = TFile::Open("NuanceCrossSections.root");
  TH2D* p_NUANCE_leptonmomentum_cc = static_cast<TH2D*>(p_fin_NUANCE->Get("leptonmomentum_cc"));
  TH2D* p_NUANCE_leptonmomentum_nc = static_cast<TH2D*>(p_fin_NUANCE->Get("leptonmomentum_nc"));
  TH2D* p_NUANCE_leptoncostheta_cc = static_cast<TH2D*>(p_fin_NUANCE->Get("leptoncostheta_cc"));
  TH2D* p_NUANCE_leptoncostheta_nc = static_cast<TH2D*>(p_fin_NUANCE->Get("leptoncostheta_nc"));
  p_NUANCE_leptonmomentum_cc->SetDirectory(0);
  p_NUANCE_leptonmomentum_nc->SetDirectory(0);
  p_NUANCE_leptoncostheta_cc->SetDirectory(0);
  p_NUANCE_leptoncostheta_nc->SetDirectory(0);
  p_fin_NUANCE->Close();


  // Load GENIE Cross Sections
  TFile* p_fin_GENIE_cc = TFile::Open("GENIECrossSections_CC.root");
  TH2D* p_GENIE_leptonmomentum_cc = static_cast<TH2D*>(p_fin_GENIE_cc->Get("leptonmomentum_cc"));
  TH2D* p_GENIE_leptoncostheta_cc = static_cast<TH2D*>(p_fin_GENIE_cc->Get("leptoncostheta_cc"));
  p_GENIE_leptonmomentum_cc->SetDirectory(0);
  p_GENIE_leptoncostheta_cc->SetDirectory(0);
  p_fin_GENIE_cc->Close();

  TFile* p_fin_GENIE_nc = TFile::Open("GENIECrossSections_NC.root");
  TH2D* p_GENIE_leptonmomentum_nc = static_cast<TH2D*>(p_fin_GENIE_nc->Get("leptonmomentum_nc"));
  TH2D* p_GENIE_leptoncostheta_nc = static_cast<TH2D*>(p_fin_GENIE_nc->Get("leptoncostheta_nc"));
  p_GENIE_leptonmomentum_nc->SetDirectory(0);
  p_GENIE_leptoncostheta_nc->SetDirectory(0);
  p_fin_GENIE_nc->Close();

  p_NUANCE_leptonmomentum_cc->Divide(p_GENIE_leptonmomentum_cc);
  p_NUANCE_leptonmomentum_nc->Divide(p_GENIE_leptonmomentum_nc);
  p_NUANCE_leptoncostheta_cc->Divide(p_GENIE_leptoncostheta_cc);
  p_NUANCE_leptoncostheta_nc->Divide(p_GENIE_leptoncostheta_nc);
  
  TFile* p_ratios = new TFile("Ratios.root","RECREATE");
  p_ratios->cd();  
  p_NUANCE_leptonmomentum_cc->Write("leptonmomentum_cc"); 
  p_NUANCE_leptonmomentum_nc->Write("leptonmomentum_nc"); 
  p_NUANCE_leptoncostheta_cc->Write("leptoncostheta_cc"); 
  p_NUANCE_leptoncostheta_nc->Write("leptoncostheta_nc"); 
  p_ratios->Close();

}
