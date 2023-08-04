//////////////////////////////////////////
// 07/21/2023 Jihee Kim (jkim11@bnl.gov)
// Calculate resolutions of particles
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

void resolution(const char *inputFile)
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
  // Get total number of events
  Long64_t allEntries = treeReader->GetEntries();
  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Electron *electron;
  // Book histograms
  TH1* hElectronDeltaPT = new TH1D("hElectronDeltaPT",";(p_{T}^{rec}-p_{T}^{gen})/p_{T}^{gen};Number of electrons",100, -0.1, 0.1);
  TH1* hElectronDeltaEta = new TH1D("hElectronDeltaEta",";(#eta^{rec}-#eta^{gen})/#eta^{gen};Number of electrons", 100, -0.1, 0.1);

  Long64_t entry;
  Int_t i;

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
  
  cout << "** Done..." << endl;

  delete treeReader;
  delete chain;
}
