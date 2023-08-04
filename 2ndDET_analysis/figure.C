//////////////////////////////////////////
// 08/03/2023 Jihee Kim (jkim11@bnl.gov)
// Draw figures
//////////////////////////////////////////

void figure()
{
  // Setting for figures
  TStyle* kStyle = new TStyle("kStyle","Kim's Style");
  kStyle->SetOptStat(0);
  kStyle->SetOptTitle(0);
  kStyle->SetOptFit(1);
  kStyle->SetStatColor(0);
  kStyle->SetStatW(0.15);
  kStyle->SetStatH(0.10);
  kStyle->SetStatX(0.85);
  kStyle->SetStatY(0.9);
  kStyle->SetStatBorderSize(1);
  kStyle->SetLabelSize(0.045,"xyz");
  kStyle->SetTitleSize(0.050,"xyz");
  kStyle->SetTitleOffset(1.2,"x");
  kStyle->SetTitleOffset(1.3,"y");
  kStyle->SetTitleOffset(1.2,"z");
  kStyle->SetLineWidth(2);
  kStyle->SetTitleFont(42,"xyz");
  kStyle->SetLabelFont(42,"xyz");
  kStyle->SetCanvasDefW(500);
  kStyle->SetCanvasDefH(500);
  kStyle->SetCanvasColor(0);
  kStyle->SetPadTickX(1);
  kStyle->SetPadTickY(1);
  kStyle->SetPadGridX(1);
  kStyle->SetPadGridY(1);
  kStyle->SetPadLeftMargin(0.15);
  kStyle->SetPadRightMargin(0.15);
  kStyle->SetPadTopMargin(0.1);
  kStyle->SetPadBottomMargin(0.15);
  TGaxis::SetMaxDigits(3);
  gStyle->SetPalette(1);
  gROOT->SetStyle("kStyle");

  // Open a ROOT file and get a histogram.
  TFile *f1 = new TFile("./plots/figures_CC_DIS_ATHENA_2T_1M_20230802_output.root");
  TFile *f2 = new TFile("./plots/figures_CC_DIS_ATHENA_3T_1M_20230802_output.root");
  
  TH1D *hTrackPAcc2T = f1->Get<TH1D>("hTrackPAcc");
  TH1D *hTrackPAcc3T = f2->Get<TH1D>("hTrackPAcc");

  TH1D *hTrackPtAcc2T = f1->Get<TH1D>("hTrackPtAcc");
  TH1D *hTrackPtAcc3T = f2->Get<TH1D>("hTrackPtAcc");

  // Plot figures
  TCanvas *cnv1 = new TCanvas("cnv1", "cnv1");
  hTrackPAcc3T->SetLineColor(kBlack); 
  hTrackPAcc2T->SetLineColor(kRed); 
  hTrackPAcc3T->Draw();
  hTrackPAcc2T->Draw("SAME");
  TLegend* leg1 = new TLegend(0.2,0.7,0.4,0.85);
  leg1->SetBorderSize(1);
  leg1->AddEntry(hTrackPAcc3T,"3T","l");
  leg1->AddEntry(hTrackPAcc2T,"2T","l");
  leg1->Draw();    
  cnv1->SaveAs("./plots/hTrackPAccDiffMag.png"); 

  TCanvas *cnv2 = new TCanvas("cnv2", "cnv2");
  hTrackPtAcc3T->SetLineColor(kBlack);
  hTrackPtAcc2T->SetLineColor(kRed);
  hTrackPtAcc3T->Draw();
  hTrackPtAcc2T->Draw("SAME");
  TLegend* leg2 = new TLegend(0.2,0.7,0.4,0.85);
  leg2->SetBorderSize(1);
  leg2->AddEntry(hTrackPtAcc3T,"3T","l");
  leg2->AddEntry(hTrackPtAcc2T,"2T","l");
  leg2->Draw();    
  cnv2->SaveAs("./plots/hTrackPtAccDiffMag.png"); 

  cout << "** Done..." << endl;
}
