/*
This macro shows how to access the particle-level reference for reconstructed objects.
It is also shown how to loop over the jet constituents.

root -l examples/Example3.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

//------------------------------------------------------------------------------

struct TestPlots
{
  TH1 *fhistDen;
  TH1 *fhistNum;

};

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  TLegend *legend;
  TPaveText *comment;

  plots->fhistDen = result->AddHist1D(
    "den", "",
    "", "",
    200, -5.0, 0.0);

  plots->fhistNum = result->AddHist1D(
    "num", "",
    "", "",
    200, -5.0, 0.0);
}

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchJet = treeReader->UseBranch("FatJet");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle, *gentrack;
  Photon *photon;

  Track *track;
  Tower *tower;

  Jet *jet;
  TObject *object;

  TLorentzVector momentum;

  Long64_t entry;

  Int_t i, j, k, pdgCode;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    if(entry %100 == 0) cout<<entry<< endl;

    jet = (Jet*) branchJet->At(0);

    // Loop over all genparticles
    for(k = 0; k < branchParticle->GetEntriesFast(); ++k)
    {

      particle = (GenParticle*) branchParticle->At(k);
      if (particle->Status != 1)  continue;
      if (particle->Charge == 0)  continue;
      if (particle->PT < 1.0)  continue;

      plots->fhistDen->Fill(TMath::Log10(particle->P4().DeltaR(jet->P4())));

      // Loop over all jet's constituents
      for(j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
      {
        object = jet->Constituents.At(j);

        // Check if the constituent is accessible
        if(object == 0) continue;

        if(object->IsA() == Track::Class())
        {
          track = (Track*) object;
          gentrack = (GenParticle*) track->Particle.GetObject();
          if(gentrack->GetUniqueID() == particle->GetUniqueID())
          {
            //cout<<"found: "<<gentrack->GetUniqueID()<<endl;
            plots->fhistNum->Fill(TMath::Log10(particle->P4().DeltaR(jet->P4())));
          }
        }
      }
    }
  }
  //plots->fhistNum->Divide(plots->fhistDen);
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
}

//------------------------------------------------------------------------------

void DenseTracks(const char *inputFile)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  TestPlots *plots = new TestPlots;

  BookHistograms(result, plots);

  AnalyseEvents(treeReader, plots);

  PrintHistograms(result, plots);

  result->Write("results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
