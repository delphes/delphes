//////////////////////////////////////////
// 08/03/2023 Jihee Kim (jkim11@bnl.gov)
// Draw tracking definition
//////////////////////////////////////////

double CommonTrackingEfficiency(double eta, double pt);
double CommonTrackingResolution(double eta, double pt);

void tracking_def()
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

  TH2* hTrackingEfficiency1GeV = new TH2D("hTrackingEfficiency1GeV", ";#eta; Tracking efficiency",100,-5.,5.,100,0.,1.2);
  TH2* hTrackingEfficiency5GeV = new TH2D("hTrackingEfficiency5GeV", ";#eta; Tracking efficiency",100,-5.,5.,100,0.,1.2);
  TH2* hTrackingEfficiency10GeV = new TH2D("hTrackingEfficiency10GeV", ";#eta; Tracking efficiency",100,-5.,5.,100,0.,1.2);

  TH2* hTrackingResolutionEtaN3p3 = new TH2D("hTrackingResolutionEtaN3p3", ";Pt [GeV]; Tracking resolution",100,0.,100.,100,0.,1.);
  TH2* hTrackingResolutionEtaN2p3 = new TH2D("hTrackingResolutionEtaN2p3", ";Pt [GeV]; Tracking resolution",100,0.,100.,100,0.,1.);
  TH2* hTrackingResolutionEtaP1p3 = new TH2D("hTrackingResolutionEtaN1p3", ";Pt [GeV]; Tracking resolution",100,0.,100.,100,0.,1.);
  
  for (int j : {1, 5, 10})
  {
    for (int i = -40; i <= 40; ++i)
    {
      const double a = i / 10.0; // a = -4.0 to 4.0 step 0.1
      if (j == 1)
        hTrackingEfficiency1GeV->Fill(a,CommonTrackingEfficiency(a,j));
      else if (j == 5)
        hTrackingEfficiency5GeV->Fill(a,CommonTrackingEfficiency(a,j));
      else  
        hTrackingEfficiency10GeV->Fill(a,CommonTrackingEfficiency(a,j));
    }
  }

  for (int k : {1.3})
  {
    for (int l = 0; l < 100; ++l)
    {
      if (k == -3.3)
        hTrackingResolutionEtaN3p3->Fill(l,CommonTrackingResolution(k,l));
      else if (k == -2.3)
        hTrackingResolutionEtaN2p3->Fill(l,CommonTrackingResolution(k,l));
      else
        hTrackingResolutionEtaP1p3->Fill(l,CommonTrackingResolution(k,l));
    }
  }

  // Plot figures
  TCanvas *cnv1 = new TCanvas("cnv1", "cnv1");
  hTrackingEfficiency1GeV->GetXaxis()->CenterTitle(true);
  hTrackingEfficiency1GeV->GetYaxis()->CenterTitle(true);
  hTrackingEfficiency1GeV->SetMarkerStyle(7);
  hTrackingEfficiency1GeV->SetMarkerColor(kBlue);
  hTrackingEfficiency1GeV->Draw();
  cnv1->SaveAs("./plots/hTrackingEfficiency1GeV.png");

  TCanvas *cnv2 = new TCanvas("cnv2", "cnv2");
  hTrackingEfficiency5GeV->GetXaxis()->CenterTitle(true);
  hTrackingEfficiency5GeV->GetYaxis()->CenterTitle(true);
  hTrackingEfficiency5GeV->SetMarkerStyle(7);
  hTrackingEfficiency5GeV->SetMarkerColor(kBlue);
  hTrackingEfficiency5GeV->Draw();
  cnv2->SaveAs("./plots/hTrackingEfficiency5GeV.png");

  TCanvas *cnv3 = new TCanvas("cnv3", "cnv3");
  hTrackingEfficiency10GeV->GetXaxis()->CenterTitle(true);
  hTrackingEfficiency10GeV->GetYaxis()->CenterTitle(true);
  hTrackingEfficiency10GeV->SetMarkerStyle(7);
  hTrackingEfficiency10GeV->SetMarkerColor(kBlue);
  hTrackingEfficiency10GeV->Draw();
  cnv3->SaveAs("./plots/hTrackingEfficiency10GeV.png");

  TCanvas *cnv4 = new TCanvas("cnv4", "cnv4");
  hTrackingResolutionEtaN3p3->GetXaxis()->CenterTitle(true);
  hTrackingResolutionEtaN3p3->GetYaxis()->CenterTitle(true);
  hTrackingResolutionEtaN3p3->SetMarkerStyle(7);
  hTrackingResolutionEtaN3p3->SetMarkerColor(kBlue);
  hTrackingResolutionEtaN2p3->SetMarkerStyle(7);
  hTrackingResolutionEtaN2p3->SetMarkerColor(kRed);
  hTrackingResolutionEtaP1p3->SetMarkerStyle(7);
  hTrackingResolutionEtaP1p3->SetMarkerColor(kGreen);
  hTrackingResolutionEtaN3p3->Draw();
  hTrackingResolutionEtaN2p3->Draw("SAME");
  hTrackingResolutionEtaP1p3->Draw("SAME");
  cnv4->SaveAs("./plots/hTrackingResolution.png");

  cout << "** Done..." << endl;
}

double CommonTrackingEfficiency(double eta, double pt)
{
  double efficiency;
  if (eta > -3.5 && eta <= -3.0 && pt*cosh(eta) > 1.25 && pt*cosh(eta)< 6.0)
    efficiency = 0.875;
  else if (eta > -3.5 && eta <= -3.0 && pt*cosh(eta) > 6.0)
    efficiency = 0.95;
  else if (eta > -3.0 && eta <= -2.5 && pt*cosh(eta) > 0.55 && pt*cosh(eta) <2.0)
    efficiency = 0.875;
  else if (eta > -3.0 && eta <= -2.5 && pt*cosh(eta) > 2.0)
    efficiency = 0.95;
  else if (eta > -2.5 && eta <= -2.0 && pt*cosh(eta) > 0.45 && pt*cosh(eta)< 0.6)
    efficiency = 0.875;
  else if (eta > -2.5 && eta <= -2.0 && pt*cosh(eta) > 0.6)
    efficiency = 0.95;
  else if (eta > -2.0 && eta <= -1.5 && pt*cosh(eta)> 0.250 && pt*cosh(eta) < 0.500)
    efficiency = 0.875;
  else if (eta > -2.0 && eta <= -1.5 && pt*cosh(eta) > 0.500)
    efficiency = 0.95;
  else if (eta > -1.5 && eta <= -1.0 && pt*cosh(eta) > 0.150 && pt*cosh(eta) < 0.300)
    efficiency = 0.86;
  else if (eta > -1.5 && eta <= -1.0 && pt*cosh(eta) > 0.300)
    efficiency = 0.92;
  else if (eta > -1.0 && eta <= -0.5 && pt*cosh(eta) > 0.150 && pt*cosh(eta) < 0.200)
    efficiency = 0.89;
  else if (eta > -1.0 && eta <= -0.5 && pt*cosh(eta) > 0.200)
    efficiency = 0.98;
  else if (eta > -0.5 && eta <= 0.0 && pt*cosh(eta) > 0.150 && pt*cosh(eta) < 0.200)
    efficiency = 0.89;
  else if (eta > -0.5 && eta <= 0.0 && pt*cosh(eta) > 0.200)
    efficiency = 0.98;
  else if (eta > 0.0 && eta <= 0.5 && pt*cosh(eta) > 0.150 && pt*cosh(eta) < 0.200)
    efficiency = 0.89;
  else if (eta > 0.0 && eta <= 0.5 && pt*cosh(eta) > 0.200)
    efficiency = 0.98;
  else if (eta > 0.5 && eta <= 1.0 && pt*cosh(eta) > 0.150 && pt*cosh(eta) < 0.200)
    efficiency = 0.89;
  else if (eta > 0.5 && eta <= 1.0 && pt*cosh(eta) > 0.200)
    efficiency = 0.98;
  else if (eta > 1.0 && eta <= 1.5 && pt*cosh(eta) > 0.150 && pt*cosh(eta) < 0.200)
    efficiency = 0.86;
  else if (eta > 1.0 && eta <= 1.5 && pt*cosh(eta) > 0.200)
    efficiency = 0.92;
  else if (eta > 1.5 && eta <= 2.0 && pt*cosh(eta) > 0.250 && pt*cosh(eta) < 0.500)
    efficiency = 0.89;
  else if (eta > 1.5 && eta <= 2.0 && pt*cosh(eta) > 0.500)
    efficiency = 0.98;
  else if (eta > 2.0 && eta <= 2.5 && pt*cosh(eta) > 0.350 && pt*cosh(eta) < 0.700)
    efficiency = 0.88;
  else if (eta > 2.0 && eta <= 2.5 && pt*cosh(eta) > 0.700)
    efficiency = 0.97;
  else if (eta > 2.5 && eta <= 3.0 && pt*cosh(eta) > 0.550 && pt*cosh(eta) < 2.0)
    efficiency = 0.87;
  else if (eta > 2.5 && eta <= 3.0 && pt*cosh(eta) > 2.0)
    efficiency = 0.95;
  else if (eta > 3.0 && eta <= 3.5 && pt*cosh(eta) > 0.850 && pt*cosh(eta) < 4.0)
    efficiency = 0.87;
  else if (eta > 3.0 && eta <= 3.5 && pt*cosh(eta) > 4.0)
    efficiency = 0.95;
  else if (abs(eta) > 3.5)
    efficiency = 0.0;
  else
    efficiency = 0.0;

  return efficiency;
}

double CommonTrackingResolution(double eta, double pt)
{
  double resolution;
  if (eta<=-3.0 && eta>-3.5)
    resolution = sqrt( pow((1.841e-2),2) + pow((pt*cosh(eta)*7.1e-4),2)  );
  else if (eta<=-2.5 && eta>-3.0)
    resolution = sqrt( pow((1.080e-2),2) + pow((pt*cosh(eta)*1.7e-4),2)  );
  else if (eta<=-2.0 && eta>-2.5)
    resolution = sqrt( pow((6.33e-3),2) + pow((pt*cosh(eta)*0.0),2)  );
  else if (eta<=-1.5 && eta>-2.0)
    resolution = sqrt( pow((4.76e-3),2) + pow((pt*cosh(eta)*1.1e-4),2)  );
  else if (eta<=-1.0 && eta>-1.5)
    resolution = sqrt( pow((4.33e-3),2) + pow((pt*cosh(eta)*1.6e-4),2)  );
  else if (eta<=-0.5 && eta>-1.0)
    resolution = sqrt( pow((3.98e-3),2) + pow((pt*cosh(eta)*5.0e-4),2)  );
  else if (eta<= 0.0 && eta>-0.5)
    resolution = sqrt( pow((3.53e-3),2) + pow((pt*cosh(eta)*5.9e-4),2)  );
  else if (eta<=0.5 && eta>0)
    resolution = sqrt( pow((3.50e-3),2) + pow((pt*cosh(eta)*5.9e-4),2)  );
  else if (eta<=1.0 && eta>0.5)
    resolution = sqrt( pow((4.01e-3),2) + pow((pt*cosh(eta)*5.0e-4),2)   );
  else if (eta<=1.5 && eta>1.0)
    resolution = sqrt( pow((4.14e-3),2) + pow((pt*cosh(eta)*1.5e-4),2)   );
  else if (eta<=2.0 && eta>1.5)
    resolution = sqrt( pow((4.66e-3),2) + pow((pt*cosh(eta)*1.1e-4),2)   );
  else if (eta<=2.5 && eta>2.0)
    resolution = sqrt( pow((6.38e-3),2) + pow((pt*cosh(eta)*1.3e-4),2)   );
  else if (eta<=3.0 && eta>2.5)
    resolution = sqrt( pow((1.089e-2),2) + pow((pt*cosh(eta)*1.1e-4),2)   );
  else if (eta<=3.5 && eta>3.0)
    resolution = sqrt( pow((1.905e-2),2) + pow((pt*cosh(eta)*3.1e-4),2)  ); 
  else
    resolution = 0.0;
  
  return resolution;
}
