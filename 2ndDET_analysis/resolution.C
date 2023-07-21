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

//------------------------------------------------------------------------------

void resolution(const char *inputFile)
{
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

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Electron *electron;

  TH1* hElectronDeltaPT = new TH1D("hElectronDeltaPT",";(p_{T}^{rec}-p_{T}^{gen})/p_{T}^{gen};Number of electrons",100, -0.1, 0.1);
  TH1* hElectronDeltaEta = new TH1D("hElectronDeltaEta",";(#eta^{rec}-#eta^{gen})/#eta^{gen};Number of electrons", 100, -0.1, 0.1);

  Long64_t entry;

  Int_t i;


  for(entry = 0; entry < allEntries; ++entry)
  {
    treeReader->ReadEntry(entry);
    for(i = 0; i < branchElectron->GetEntriesFast(); ++i)
      {
        electron = (Electron*) branchElectron->At(i);
        particle = (GenParticle*) electron->Particle.GetObject();

        hElectronDeltaPT->Fill((particle->PT - electron->PT)/particle->PT);
        hElectronDeltaEta->Fill((particle->Eta - electron->Eta)/particle->Eta);
      }
  }

  TCanvas *cnv1 = new TCanvas("cnv1", "cnv1");
  hElectronDeltaPT->GetXaxis()->CenterTitle(true);
  hElectronDeltaPT->GetYaxis()->CenterTitle(true);
  hElectronDeltaPT->SetLineWidth(2);
  //double up_fit   = hElectronDeltaPT->GetMean() + 1.0*hElectronDeltaPT->GetStdDev();
  //double down_fit = hElectronDeltaPT->GetMean() - 1.0*hElectronDeltaPT->GetStdDev();
  //hElectronDeltaPT->Fit("gaus", "", "", down_fit, up_fit);
  hElectronDeltaPT->Fit("gaus");
  TF1 *gaus = hElectronDeltaPT->GetFunction("gaus");
  gaus->SetLineWidth(2);
  gaus->SetLineColor(kRed);
  hElectronDeltaPT->Draw();
  cnv1->SaveAs("./2ndDET_analysis/hElectronDeltaPT.png");

  cout << "** Done..." << endl;

  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
