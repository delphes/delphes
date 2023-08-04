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

void track(const char *inputFile, const char *outputFile)
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
  TH1* hTrackPGen = new TH1D("hTrackPGen", ";Generated P [GeV];Number of tracks",100,0.,200.);
  TH1* hTrackPRec = new TH1D("hTrackPRec", ";Reconstructed P [GeV];Number of tracks",100,0.,200.);
  TH1* hTrackPAcc = new TH1D("hTrackPAcc", ";P [GeV];P Acceptance = Rec/Gen",100,0.,200.);
  TH1* hTrackDeltaP = new TH1D("hTrackDeltaP",";(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 

  TH1* hTrackPtGen = new TH1D("hTrackPtGen", ";Generated Pt [GeV];Number of tracks",100,0.,60.);
  TH1* hTrackPtRec = new TH1D("hTrackPtRec", ";Reconstructed Pt [GeV];Number of tracks",100,0.,60.);
  TH1* hTrackPtAcc = new TH1D("hTrackPtAcc", ";Pt [GeV];Pt Acceptance = Rec/Gen",100,0.,60.);
  
  TH1* hTrackEtaGen = new TH1D("hTrackEtaGen", ";Generated #eta;Number of tracks",100,-5.,5.);
  TH1* hTrackEtaRec = new TH1D("hTrackEtaRec", ";Reconstructed #eta;Number of tracks",100,-5.,5.);
  TH1* hTrackEtaAcc = new TH1D("hTrackEtaAcc", ";#eta;#eta Acceptance = Rec/Gen",100,-5.,5.);
  TH1* hTrackDeltaEta = new TH1D("hTrackDeltaEta",";(#eta^{rec}-#eta^{gen})/#eta^{gen};Number of tracks", 100, -0.1, 0.1);

  TH2* hTrackPvsEtaGen = new TH2D("hTrackPvsEtaGen", ";#eta^{gen};Generated P [GeV]",100,-5.,5.,100,0.,200.);
  TH2* hTrackPvsEtaRec = new TH2D("hTrackPvsEtaRec", ";#eta^{rec};Reconstructed P [GeV]",100,-5.,5.,100,0.,200.);

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
	
        // Track momentum
        hTrackPGen->Fill(particle->P);
        hTrackPRec->Fill(track->P);
        hTrackDeltaP->Fill((track->P - particle->P)/particle->P);
        // Track transverse momentum
        hTrackPtGen->Fill(particle->PT);
        hTrackPtRec->Fill(track->PT);
        // Track pseudorapidity
        hTrackEtaGen->Fill(particle->Eta);
        hTrackEtaRec->Fill(track->Eta);
        hTrackDeltaEta->Fill((track->Eta - particle->Eta)/particle->Eta);
        // Track momentum versus pseudorapidity
        hTrackPvsEtaGen->Fill(particle->Eta,particle->P);
        hTrackPvsEtaRec->Fill(track->Eta,track->P);
      }
    }
  }

  // Output file to save all figures
  TFile* output = new TFile(outputFile,"RECREATE");
  
  // Plot figures
  TCanvas *cnv1 = new TCanvas("cnv1", "cnv1");
  cnv1->SetLogy(1);
  hTrackPGen->GetXaxis()->CenterTitle(true);
  hTrackPGen->GetYaxis()->CenterTitle(true);
  hTrackPGen->SetLineWidth(2);
  hTrackPGen->Draw();
  cnv1->SaveAs("./plots/hTrackPGen.png");

  TCanvas *cnv2 = new TCanvas("cnv2", "cnv2");
  cnv2->SetLogy(1);
  hTrackPRec->GetXaxis()->CenterTitle(true);
  hTrackPRec->GetYaxis()->CenterTitle(true);
  hTrackPRec->SetLineWidth(2);
  hTrackPRec->Draw();
  cnv2->SaveAs("./plots/hTrackPRec.png");

  TCanvas *cnv3 = new TCanvas("cnv3", "cnv3");
  Int_t NbP = hTrackPGen->GetNbinsX();
  for (Int_t k = 1; k < NbP + 1; k++)
  {
    Float_t pgen = hTrackPGen->GetBinContent(k);
    Float_t prec = hTrackPRec->GetBinContent(k);
    Float_t PAcc = 0.0;
    if (pgen > 0.0) PAcc = prec / pgen;
    hTrackPAcc->SetBinContent(k, PAcc);
  }
  hTrackPAcc->GetYaxis()->SetRangeUser(0.,3.2); // CC_DIS
  //hTrackPAcc->GetYaxis()->SetRangeUser(0.,4.2); // NC_DIS
  hTrackPAcc->GetXaxis()->CenterTitle(true);
  hTrackPAcc->GetYaxis()->CenterTitle(true);
  hTrackPAcc->SetLineWidth(2);
  hTrackPAcc->Draw(); 
  hTrackPAcc->Write(); 
  cnv3->SaveAs("./plots/hTrackPAcc.png");

  TCanvas *cnv4 = new TCanvas("cnv4", "cnv4");
  hTrackDeltaP->GetXaxis()->CenterTitle(true);
  hTrackDeltaP->GetYaxis()->CenterTitle(true);
  hTrackDeltaP->SetLineWidth(2);
  hTrackDeltaP->Fit("gaus");
  TF1 *gaus_hTrackDeltaP = hTrackDeltaP->GetFunction("gaus");
  gaus_hTrackDeltaP->SetLineWidth(2);
  gaus_hTrackDeltaP->SetLineColor(kRed);
  hTrackDeltaP->Draw();
  cnv4->SaveAs("./plots/hTrackDeltaP.png");

  TCanvas *cnv5 = new TCanvas("cnv5", "cnv5");
  cnv5->SetLogy(1);
  hTrackPtGen->GetXaxis()->CenterTitle(true);
  hTrackPtGen->GetYaxis()->CenterTitle(true);
  hTrackPtGen->SetLineWidth(2);
  hTrackPtGen->Draw();
  cnv5->SaveAs("./plots/hTrackPtGen.png");

  TCanvas *cnv6 = new TCanvas("cnv6", "cnv6");
  cnv6->SetLogy(1);
  hTrackPtRec->GetXaxis()->CenterTitle(true);
  hTrackPtRec->GetYaxis()->CenterTitle(true);
  hTrackPtRec->SetLineWidth(2);
  hTrackPtRec->Draw();
  cnv6->SaveAs("./plots/hTrackPtRec.png");
  
  TCanvas *cnv7 = new TCanvas("cnv7", "cnv7");
  Int_t NbPt = hTrackPtGen->GetNbinsX();
  for (Int_t j = 1; j < NbPt + 1; j++)
  {
    Float_t ptgen = hTrackPtGen->GetBinContent(j);
    Float_t ptrec = hTrackPtRec->GetBinContent(j);
    Float_t PtAcc = 0.0;
    if (ptgen > 0.0) PtAcc = ptrec / ptgen;
    hTrackPtAcc->SetBinContent(j, PtAcc);
  }
  hTrackPAcc->GetYaxis()->SetRangeUser(0.,5.2); // CC_DIS
  //hTrackPAcc->GetYaxis()->SetRangeUser(0.,2.8); // NC_DIS
  hTrackPtAcc->GetXaxis()->CenterTitle(true);
  hTrackPtAcc->GetYaxis()->CenterTitle(true);
  hTrackPtAcc->SetLineWidth(2);
  hTrackPtAcc->Draw(); 
  hTrackPtAcc->Write(); 
  cnv7->SaveAs("./plots/hTrackPtAcc.png");

  TCanvas *cnv8 = new TCanvas("cnv8", "cnv8");
  hTrackEtaGen->GetXaxis()->CenterTitle(true);
  hTrackEtaGen->GetYaxis()->CenterTitle(true);
  hTrackEtaGen->SetLineWidth(2);
  hTrackEtaGen->Draw();
  cnv8->SaveAs("./plots/hTrackEtaGen.png");

  TCanvas *cnv9 = new TCanvas("cnv9", "cnv9");
  hTrackEtaRec->GetXaxis()->CenterTitle(true);
  hTrackEtaRec->GetYaxis()->CenterTitle(true);
  hTrackEtaRec->SetLineWidth(2);
  hTrackEtaRec->Draw();
  cnv9->SaveAs("./plots/hTrackEtaRec.png");

  TCanvas *cnv10 = new TCanvas("cnv10", "cnv10");
  Int_t NbEta = hTrackEtaGen->GetNbinsX();
  for (Int_t l = 1; l < NbEta + 1; l++)
  {
    Float_t etagen = hTrackEtaGen->GetBinContent(l);
    Float_t etarec = hTrackEtaRec->GetBinContent(l);
    Float_t EtaAcc = 0.0;
    if (etagen > 0.0) EtaAcc = etarec / etagen;
    hTrackEtaAcc->SetBinContent(l, EtaAcc);
  }
  hTrackEtaAcc->GetXaxis()->CenterTitle(true);
  hTrackEtaAcc->GetYaxis()->CenterTitle(true);
  hTrackEtaAcc->SetLineWidth(2);
  hTrackEtaAcc->Draw(); 
  cnv10->SaveAs("./plots/hTrackEtaAcc.png");

  TCanvas *cnv11 = new TCanvas("cnv11", "cnv11");
  hTrackDeltaEta->GetXaxis()->CenterTitle(true);
  hTrackDeltaEta->GetYaxis()->CenterTitle(true);
  hTrackDeltaEta->SetLineWidth(2);
  hTrackDeltaEta->Fit("gaus");
  TF1 *gaus_hTrackDeltaEta = hTrackDeltaEta->GetFunction("gaus");
  gaus_hTrackDeltaEta->SetLineWidth(2);
  gaus_hTrackDeltaEta->SetLineColor(kRed);
  hTrackDeltaEta->Draw();
  cnv11->SaveAs("./plots/hTrackDeltaEta.png");

  TCanvas *cnv12 = new TCanvas("cnv12", "cnv12");
  hTrackPvsEtaGen->GetXaxis()->CenterTitle(true);
  hTrackPvsEtaGen->GetYaxis()->CenterTitle(true);
  hTrackPvsEtaGen->SetLineWidth(2);
  hTrackPvsEtaGen->Draw("COLZ");
  cnv12->SaveAs("./plots/hTrackPvsEtaGen.png");

  TCanvas *cnv13 = new TCanvas("cnv13", "cnv13");
  hTrackPvsEtaRec->GetXaxis()->CenterTitle(true);
  hTrackPvsEtaRec->GetYaxis()->CenterTitle(true);
  hTrackPvsEtaRec->SetLineWidth(2);
  hTrackPvsEtaRec->Draw("COLZ");
  cnv13->SaveAs("./plots/hTrackPvsEtaRec.png");

  cout << "** Done..." << endl;

  delete treeReader;
  delete chain;
}
