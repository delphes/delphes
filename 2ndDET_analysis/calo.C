//////////////////////////////////////////
// 08/10/2023 Jihee Kim (jkim11@bnl.gov)
//  
//////////////////////////////////////////

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes.so)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

void calo(const char *inputFile)
{
  // Setting for figures
  TStyle* kStyle = new TStyle("kStyle","Kim's Style");
  kStyle->SetOptStat("emr");
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

  gSystem->Load("libDelphes");
  // Create chain of root trees
  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  // Get pointers to branches used in this analysis
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");
  // Get total number of events
  Long64_t allEntries = treeReader->GetEntries();
  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Electron *electron;
  Photon *photon;
  Tower *tower;
  // Book histograms
  TH1* hElectronDeltaPT = new TH1D("hElectronDeltaPT",";(p_{T}^{rec}-p_{T}^{gen})/p_{T}^{gen};Number of electrons",100, -0.1, 0.1);
  TH1* hElectronDeltaEta = new TH1D("hElectronDeltaEta",";(#eta^{rec}-#eta^{gen})/#eta^{gen};Number of electrons", 100, -0.1, 0.1);
  
  TH1* hPhotonDeltaPT = new TH1D("hPhotonDeltaPT",";(p_{T}^{rec}-p_{T}^{gen})/p_{T}^{gen};Number of photons",100, -0.5, 0.5);
  TH1* hPhotonDeltaEta = new TH1D("hPhotonDeltaEta",";(#eta^{rec}-#eta^{gen})/#eta^{gen};Number of photons", 100, -0.3, 0.3);
  TH1* hPhotonDeltaE = new TH1D("hPhotonDeltaE",";(E^{rec}-E^{gen})/E^{gen};Number of photons", 100, -1.0, 1.0);

  TH2* hTowerEtaPhi = new TH2D("hTowerEtaPhi",";#eta; #phi [rad]",100,-5.,5.,100,-3.14,3.14);
  TH2* hTowerEtaE = new TH2D("hTowerEtaE",";#eta; E [GeV]",100,-5.,5.,100,0.,200.);

  Long64_t entry;
  Int_t i, j, k;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    // Loop over all electrons in event
    for(i = 0; i < branchElectron->GetEntriesFast(); ++i)
    {
      electron = (Electron*) branchElectron->At(i);
      particle = (GenParticle*) electron->Particle.GetObject();

      hElectronDeltaPT->Fill((electron->PT - particle->PT)/particle->PT);
      hElectronDeltaEta->Fill((electron->Eta - particle->Eta)/particle->Eta);
    }
    // Loop over all photons in event
    for(j = 0; j < branchPhoton->GetEntriesFast(); ++j)
    {
      photon = (Photon*) branchPhoton->At(j);
      
      // Skip photons with references to multiple particles
      if(photon->Particles.GetEntriesFast() != 1) continue;
      particle = (GenParticle*) photon->Particles.At(0);
      
      hPhotonDeltaPT->Fill((photon->PT - particle->PT)/particle->PT);
      hPhotonDeltaEta->Fill((photon->Eta - particle->Eta)/particle->Eta);
      hPhotonDeltaE->Fill((photon->E - particle->E)/particle->E);                                         
    }
    // Loop over all tower in event
    for(k = 0; k < branchTower->GetEntriesFast(); ++k)
    {
      tower = (Tower*) branchTower->At(k);
      hTowerEtaPhi->Fill(tower->Eta,tower->Phi);
      hTowerEtaE->Fill(tower->Eta,tower->E);
    }
  }

  // Plot figures
  TCanvas *cnv1 = new TCanvas("cnv1", "cnv1");
  hElectronDeltaPT->GetXaxis()->CenterTitle(true);
  hElectronDeltaPT->GetYaxis()->CenterTitle(true);
  hElectronDeltaPT->SetLineWidth(2);
  hElectronDeltaPT->Fit("gaus");
  TF1 *gaus_hElectronDeltaPT = hElectronDeltaPT->GetFunction("gaus");
  gaus_hElectronDeltaPT->SetLineWidth(2);
  gaus_hElectronDeltaPT->SetLineColor(kRed);
  hElectronDeltaPT->Draw();
  cnv1->SaveAs("./plots/hElectronDeltaPT.png");

  TCanvas *cnv2 = new TCanvas("cnv2", "cnv2");
  hElectronDeltaEta->GetXaxis()->CenterTitle(true);
  hElectronDeltaEta->GetYaxis()->CenterTitle(true);
  hElectronDeltaEta->SetLineWidth(2);
  hElectronDeltaEta->Fit("gaus");
  TF1 *gaus_hElectronDeltaEta = hElectronDeltaEta->GetFunction("gaus");
  gaus_hElectronDeltaEta->SetLineWidth(2);
  gaus_hElectronDeltaEta->SetLineColor(kRed);
  hElectronDeltaEta->Draw();
  cnv2->SaveAs("./plots/hElectronDeltaEta.png");
 
  TCanvas *cnv3 = new TCanvas("cnv3", "cnv3");
  hPhotonDeltaPT->GetXaxis()->CenterTitle(true);
  hPhotonDeltaPT->GetYaxis()->CenterTitle(true);
  hPhotonDeltaPT->SetLineWidth(2);
  hPhotonDeltaPT->Fit("gaus");
  TF1 *gaus_hPhotonDeltaPT = hPhotonDeltaPT->GetFunction("gaus");
  gaus_hPhotonDeltaPT->SetLineWidth(2);
  gaus_hPhotonDeltaPT->SetLineColor(kRed);
  hPhotonDeltaPT->Draw();
  cnv3->SaveAs("./plots/hPhotonDeltaPT.png");

  TCanvas *cnv4 = new TCanvas("cnv4", "cnv4");
  hPhotonDeltaEta->GetXaxis()->CenterTitle(true);
  hPhotonDeltaEta->GetYaxis()->CenterTitle(true);
  hPhotonDeltaEta->SetLineWidth(2);
  hPhotonDeltaEta->Fit("gaus");
  TF1 *gaus_hPhotonDeltaEta = hPhotonDeltaEta->GetFunction("gaus");
  gaus_hPhotonDeltaEta->SetLineWidth(2);
  gaus_hPhotonDeltaEta->SetLineColor(kRed);
  hPhotonDeltaEta->Draw();
  cnv4->SaveAs("./plots/hPhotonDeltaEta.png");

  TCanvas *cnv5 = new TCanvas("cnv5", "cnv5");
  hPhotonDeltaE->GetXaxis()->CenterTitle(true);
  hPhotonDeltaE->GetYaxis()->CenterTitle(true);
  hPhotonDeltaE->SetLineWidth(2);
  hPhotonDeltaE->Fit("gaus");
  TF1 *gaus_hPhotonDeltaE = hPhotonDeltaE->GetFunction("gaus");
  gaus_hPhotonDeltaE->SetLineWidth(2);
  gaus_hPhotonDeltaE->SetLineColor(kRed);
  hPhotonDeltaE->Draw();
  cnv5->SaveAs("./plots/hPhotonDeltaE.png");

  TCanvas *cnv6 = new TCanvas("cnv6", "cnv6");
  hTowerEtaPhi->GetXaxis()->CenterTitle(true);
  hTowerEtaPhi->GetYaxis()->CenterTitle(true);
  hTowerEtaPhi->Draw("COLZ");
  cnv6->SaveAs("./plots/hTowerEtaPhi.png");

  TCanvas *cnv7 = new TCanvas("cnv7", "cnv7");
  hTowerEtaE->GetXaxis()->CenterTitle(true);
  hTowerEtaE->GetYaxis()->CenterTitle(true);
  hTowerEtaE->Draw("COLZ");
  cnv7->SaveAs("./plots/hTowerEtaE.png");

  cout << "** Done..." << endl;

  delete treeReader;
  delete chain;
}
