//////////////////////////////////////////
// 08/02/2023 Jihee Kim (jkim11@bnl.gov)
// Tracks of particles
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

void track(const char *inputFile)
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
  TClonesArray* branchTrack = treeReader->UseBranch("Track");
  // Get total number of events
  Long64_t allEntries = treeReader->GetEntries();
  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Track *track;  
  
  // Book histograms
  TH1* hTrackPtGen = new TH1D("hTrackPtGen", ";Generated Pt [GeV];Number of tracks",100,0.,60.);
  TH1* hTrackPtRec = new TH1D("hTrackPtRec", ";Reconstructed Pt [GeV];Number of tracks",100,0.,60.);
  TH1* hTrackPtAcc = new TH1F("hTrackPtAcc", ";Pt [GeV];Pt Acceptance = Rec/Gen",100,0.,60.);

  Long64_t entry;
  Int_t i;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    //
    if (branchTrack->GetEntries() > 0)
    {
      // Loop on all tracks
      for (i = 0; i < branchTrack->GetEntries(); ++i)
      {
        track = (Track*) branchTrack->At(i);
        particle = (GenParticle*) track->Particle.GetObject();

        hTrackPtGen->Fill(particle->PT);
        hTrackPtRec->Fill(track->PT);
      }
    }
  }

  // Plot figures
  TCanvas *cnv1 = new TCanvas("cnv1", "cnv1");
  cnv1->SetLogy(1);
  hTrackPtGen->GetXaxis()->CenterTitle(true);
  hTrackPtGen->GetYaxis()->CenterTitle(true);
  hTrackPtGen->SetLineWidth(2);
  hTrackPtGen->Draw();
  cnv1->SaveAs("./plots/hTrackPtGen.png");

  TCanvas *cnv2 = new TCanvas("cnv2", "cnv2");
  cnv2->SetLogy(1);
  hTrackPtRec->GetXaxis()->CenterTitle(true);
  hTrackPtRec->GetYaxis()->CenterTitle(true);
  hTrackPtRec->SetLineWidth(2);
  hTrackPtRec->Draw();
  cnv2->SaveAs("./plots/hTrackPtRec.png");
  
  TCanvas *cnv3 = new TCanvas("cnv3", "cnv3");
  Int_t NbPt = hTrackPtGen->GetNbinsX();
  for (Int_t j = 1; j < NbPt + 1; j++)
  {
    Float_t pgen = hTrackPtGen->GetBinContent(j);
    Float_t prec = hTrackPtRec->GetBinContent(j);
    Float_t PtAcc = 0.0;
    if (pgen > 0.0) PtAcc = prec / pgen;
    hTrackPtAcc->SetBinContent(j, PtAcc);
  }
  hTrackPtAcc->GetXaxis()->CenterTitle(true);
  hTrackPtAcc->GetYaxis()->CenterTitle(true);
  hTrackPtAcc->SetLineWidth(2);
  hTrackPtAcc->Draw(); 
  cnv3->SaveAs("./plots/hTrackPtAcc.png");
  
  cout << "** Done..." << endl;

  delete treeReader;
  delete chain;
}
