/*
Simple macro showing how to access branches from the delphes output root file for Snowmass studies,
loop over events, and plot simple quantities such as the leading electron and jet pt.

root -l examples/Example7.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

//------------------------------------------------------------------------------

void Example7(const char *inputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("JetPUPPITight");
  TClonesArray *branchElectron = treeReader->UseBranch("ElectronMedium");
  TClonesArray *branchWeight = treeReader->UseBranch("Weight");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  // Book histograms
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 1000.0);
  TH1 *histElectronPT = new TH1F("ele_pt", "electron P_{T}", 100, 0.0, 1000.0);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    
    // main MC event weight
    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    Double_t w = event->Weight;
    
    // read lhe event weights
    for(Int_t i = 0; i < branchWeight->GetEntriesFast(); ++i)
    {
      Weight *weight = (Weight*) branchWeight -> At(i);
      Double_t lhe_weight = weight->Weight;
      
      //cout<<lhe_weight<<endl;
      // do stuff ...
    }

    // If event contains at least 1 jet
    if(branchJet->GetEntries() > 0)
    {
      // Take first jet
      Jet *jet = (Jet*) branchJet->At(0);

      // 0 - Loose , 1 - Medium, 2 - Tight
      Int_t wp = 1;

      Bool_t BtagOk = ( jet->BTag & (1 << wp) );
      Double_t pt = jet->PT;
      Double_t eta = TMath::Abs(jet->Eta);

      // require jet to be within acceptance and btagged with Medium WP
      if (BtagOk && pt > 30. && eta < 5.)  
        histJetPT->Fill(jet->PT, w);
    }

    // If event contains at least 1 jet
    if(branchElectron->GetEntries() > 0)
    {
      // Take first jet
      Electron *electron = (Electron*) branchElectron->At(0);

      Double_t pt = electron->PT;
      Double_t eta = TMath::Abs(electron->Eta);

      // 0.1 - Loose , 0.2 - Medium, 0.3 - Tight
      Double_t IsoCut = 0.2;
      Bool_t IsoOk = electron->IsolationVar < IsoCut;

      // Plot jet transverse momentum
      if (IsoOk && pt > 10. && eta < 5.)  
        histElectronPT->Fill(electron->PT, w);
    }
  }


  //
  TCanvas *cnv = new TCanvas("cnv", "cnv", 50, 50, 800, 500);
  cnv->Divide(2, 1);
  cnv->cd(1);
  gStyle->SetOptStat(0);

  histJetPT->Draw();

  cnv->cd(2);
  gStyle->SetOptStat(0);
  // Show resulting histograms
  histElectronPT->Draw();

  cnv->Print("example7.png", "png");

}

