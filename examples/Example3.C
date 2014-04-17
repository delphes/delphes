/*
root -l examples/Example3.C\(\"delphes_output.root\"\)
*/

//------------------------------------------------------------------------------

struct TestPlots
{
  TH1 *fElectronDeltaPT;
  TH1 *fElectronDeltaEta;

  TH1 *fPhotonDeltaPT;
  TH1 *fPhotonDeltaEta;
  TH1 *fPhotonDeltaE;

  TH1 *fMuonDeltaPT;
  TH1 *fMuonDeltaEta;

  TH1 *fTrackDeltaPT;
  TH1 *fTrackDeltaEta;

  TH1 *fJetDeltaPT;
};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  TLegend *legend;
  TPaveText *comment;

  plots->fElectronDeltaPT = result->AddHist1D(
    "electron delta pt", "(p_{T}^{particle} - p_{T}^{electron})/p_{T}^{particle}",
    "(p_{T}^{particle} - p_{T}^{electron})/p_{T}^{particle}", "number of electrons",
    100, -0.1, 0.1);

  plots->fElectronDeltaEta = result->AddHist1D(
    "electron delta eta", "(#eta^{particle} - #eta^{electron})/#eta^{particle}",
    "(#eta^{particle} - #eta^{electron})/#eta^{particle}", "number of electrons",
    100, -0.1, 0.1);

  plots->fPhotonDeltaPT = result->AddHist1D(
    "photon delta pt", "(p_{T}^{particle} - p_{T}^{photon})/p_{T}^{particle}",
    "(p_{T}^{particle} - p_{T}^{photon})/p_{T}^{particle}", "number of photons",
    100, -0.1, 0.1);

  plots->fPhotonDeltaEta = result->AddHist1D(
    "photon delta eta", "(#eta^{particle} - #eta^{photon})/#eta^{particle}",
    "(#eta^{particle} - #eta^{photon})/#eta^{particle}", "number of photons",
    100, -0.1, 0.1);

  plots->fPhotonDeltaE = result->AddHist1D(
    "photon delta energy", "(E^{particle} - E^{photon})/E^{particle}",
    "(E^{particle} - E^{photon})/E^{particle}", "number of photons",
    100, -0.1, 0.1);

  plots->fMuonDeltaPT = result->AddHist1D(
    "muon delta pt", "(p_{T}^{particle} - p_{T}^{muon})/p_{T}^{particle}",
    "(p_{T}^{particle} - p_{T}^{muon})/p_{T}^{particle}", "number of muons",
    100, -0.1, 0.1);

  plots->fMuonDeltaEta = result->AddHist1D(
    "muon delta eta", "(#eta^{particle} - #eta^{muon})/#eta^{particle}",
    "(#eta^{particle} - #eta^{muon})/#eta^{particle}", "number of muons",
    100, -0.1, 0.1);

  plots->fTrackDeltaPT = result->AddHist1D(
    "track delta pt", "(p_{T}^{particle} - p_{T}^{track})/p_{T}^{particle}",
    "(p_{T}^{particle} - p_{T}^{track})/p_{T}^{particle}", "number of tracks",
    100, -0.1, 0.1);

  plots->fTrackDeltaEta = result->AddHist1D(
    "track delta eta", "(#eta^{particle} - #eta^{track})/#eta^{particle}",
    "(#eta^{particle} - #eta^{track})/#eta^{particle}", "number of tracks",
    100, -0.1, 0.1);

  plots->fJetDeltaPT = result->AddHist1D(
    "jet delta pt", "(p_{T}^{jet} - p_{T}^{constituents})/p_{T}^{jet}",
    "(p_{T}^{jet} - p_{T}^{constituents})/p_{T}^{jet}", "number of jets",
    100, -1.0e-1, 1.0e-1);

}

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");

  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");

  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Electron *electron;
  Photon *photon;
  Muon *muon;

  Track *track;
  Tower *tower;

  Jet *jet;
  TObject *object;

  TLorentzVector momentum;

  Float_t Eem, Ehad;
  Bool_t skip;

  Long64_t entry;

  Int_t i, j, pdgCode;

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

      plots->fElectronDeltaPT->Fill((particle->PT - electron->PT)/particle->PT);
      plots->fElectronDeltaEta->Fill((particle->Eta - electron->Eta)/particle->Eta);
    }

    // Loop over all photons in event
    for(i = 0; i < branchPhoton->GetEntriesFast(); ++i)
    {
      photon = (Photon*) branchPhoton->At(i);

      // skip photons with references to multiple particles
      if(photon->Particles.GetEntriesFast() != 1) continue;

      particle = (GenParticle*) photon->Particles.At(0);

      plots->fPhotonDeltaPT->Fill((particle->PT - photon->PT)/particle->PT);
      plots->fPhotonDeltaEta->Fill((particle->Eta - photon->Eta)/particle->Eta);
      plots->fPhotonDeltaE->Fill((particle->E - photon->E)/particle->E);
    }

    // Loop over all muons in event
    for(i = 0; i < branchMuon->GetEntriesFast(); ++i)
    {
      muon = (Muon*) branchMuon->At(i);
      particle = (GenParticle*) muon->Particle.GetObject();

      plots->fMuonDeltaPT->Fill((particle->PT - muon->PT)/particle->PT);
      plots->fMuonDeltaEta->Fill((particle->Eta - muon->Eta)/particle->Eta);
    }

    // Loop over all tracks in event
    for(i = 0; i < branchTrack->GetEntriesFast(); ++i)
    {
      track = (Track*) branchTrack->At(i);
      particle = (GenParticle*) track->Particle.GetObject();

      plots->fTrackDeltaPT->Fill((particle->PT - track->PT)/particle->PT);
      plots->fTrackDeltaEta->Fill((particle->Eta - track->Eta)/particle->Eta);
    }

    cout<<"--  New event -- "<<endl;

    // Loop over all jets in event
    for(i = 0; i < branchJet->GetEntriesFast(); ++i)
    {
      jet = (Jet*) branchJet->At(i);

      momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

      cout<<"Looping over jet constituents. Jet pt: "<<jet->PT<<", eta: "<<jet->Eta<<", phi: "<<jet->Phi<<endl;      

      // Loop over all jet's constituents
      for(j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
      {
        object = jet->Constituents.At(j);

        // Check if the constituent is accessible
        if(object == 0) continue;

        if(object->IsA() == GenParticle::Class())
        {
          particle = (GenParticle*) object;
          cout << "    GenPart pt: " << particle->PT << ", eta: " << particle->Eta << ", phi: " << particle->Phi << endl;
          momentum += particle->P4();
        }
        else if(object->IsA() == Track::Class())
        {
          track = (Track*) object;
          cout << "    Track pt: " << track->PT << ", eta: " << track->Eta << ", phi: " << track->Phi << endl;
          momentum += track->P4();
        }
        else if(object->IsA() == Tower::Class())
        {
          tower = (Tower*) object;
          cout << "    Tower pt: " << tower->ET << ", eta: " << tower->Eta << ", phi: " << tower->Phi << endl;
          momentum += tower->P4();
        }
        else if(object->IsA() == Muon::Class())
        {
          muon = (Muon*) object;
          cout << "    Muon pt: " << muon->PT << ", eta: " << muon->Eta << ", phi: " << muon->Phi << endl;
          momentum += muon->P4();
        }
      }
      plots->fJetDeltaPT->Fill((jet->PT - momentum.Pt())/jet->PT);
    }
  }
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
}

//------------------------------------------------------------------------------

void Example3(const char *inputFile)
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
