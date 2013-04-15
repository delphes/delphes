#include <stdexcept>
#include <iostream>
#include <sstream>
#include <memory>

#include <map>

#include <stdlib.h>
#include <signal.h>
#include <stdio.h>

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TObjArray.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

#include "modules/Delphes.h"
#include "classes/DelphesStream.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace std;

//---------------------------------------------------------------------------

void ConvertInput(fwlite::Event &event, DelphesFactory *factory, TObjArray *allParticleOutputArray, TObjArray *stableParticleOutputArray, TObjArray *partonOutputArray)
{
  size_t i;
  fwlite::Handle< vector< reco::GenParticle > > handleParticle;
  handleParticle.getByLabel(event, "genParticles");

  Candidate *candidate;
  TDatabasePDG *pdg;
  TParticlePDG *pdgParticle;
  Int_t pdgCode;

  Int_t pid, status;
  Double_t px, py, pz, e;
  Double_t x, y, z;

  pdg = TDatabasePDG::Instance();

  for(i = 0; i < handleParticle->size(); ++i)
  {
    const reco::GenParticle &particle = (*handleParticle)[i];

     pid = particle.pdgId();
     status = particle.status();
     px = particle.px(), py = particle.py(), pz = particle.phi(), e = particle.energy();
     x = particle.vx(), y = particle.vy(), z = particle.vz();

    if(status == 1 || status == 2)
    {
      candidate = factory->NewCandidate();

      candidate->PID = pid;
      pdgCode = TMath::Abs(candidate->PID);

      candidate->Status = status;

      pdgParticle = pdg->GetParticle(pid);
      candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
      candidate->Mass = pdgParticle ? pdgParticle->Mass() : -999.9;

      candidate->Momentum.SetPxPyPzE(px, py, pz, e);

      candidate->Position.SetXYZT(x, y, z, 0.0);

      allParticleOutputArray->Add(candidate);

      if(!pdgParticle) return;

      if(status == 1)
      {
        stableParticleOutputArray->Add(candidate);
      }
      else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15)
      {
        partonOutputArray->Add(candidate);
      }
    }
  }
}

//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "DelphesCMSFWLite";
  stringstream message;
  TFile *inputFile = 0;
  TFile *outputFile = 0;
  TStopwatch eventStopWatch;
  ExRootTreeWriter *treeWriter = 0;
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  DelphesFactory *factory = 0;
  TObjArray *allParticleOutputArray = 0, *stableParticleOutputArray = 0, *partonOutputArray = 0;
  Int_t i;
  Long64_t entry, allEntries;

  if(argc < 4)
  {
    cout << " Usage: " << appName << " config_file" << " output_file" << " input_file(s)" << endl;
    cout << " config_file - configuration file in Tcl format," << endl;
    cout << " output_file - output file in ROOT format," << endl;
    cout << " input_file(s) - input file(s) in ROOT format." << endl;
    return 1;
  }

  signal(SIGINT, SignalHandler);

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  AutoLibraryLoader::enable();

  try
  {
    outputFile = TFile::Open(argv[2], "CREATE");

    if(outputFile == NULL)
    {
      message << "can't open " << argv[2] << endl;
      throw runtime_error(message.str());
    }

    treeWriter = new ExRootTreeWriter(outputFile, "Delphes");

    confReader = new ExRootConfReader;
    confReader->ReadFile(argv[1]);

    modularDelphes = new Delphes("Delphes");
    modularDelphes->SetConfReader(confReader);
    modularDelphes->SetTreeWriter(treeWriter);

    factory = modularDelphes->GetFactory();
    allParticleOutputArray = modularDelphes->ExportArray("allParticles");
    stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
    partonOutputArray = modularDelphes->ExportArray("partons");

    modularDelphes->InitTask();

    for(i = 3; i < argc && !interrupted; ++i)
    {
      cout << "** Reading " << argv[i] << endl;

      inputFile = TFile::Open(argv[i]);

      if(inputFile == NULL)
      {
        message << "can't open " << argv[i] << endl;
        throw runtime_error(message.str());
      }

      fwlite::Event event(inputFile);

      allEntries = event.size();

      if(allEntries <= 0) continue;

      ExRootProgressBar progressBar(allEntries - 1);

      // Loop over all objects
      entry = 0;
      modularDelphes->Clear();
      treeWriter->Clear();
      for(event.toBegin(); !event.atEnd() && !interrupted; ++event)
      {
        ConvertInput(event, factory, allParticleOutputArray, stableParticleOutputArray, partonOutputArray);
        modularDelphes->ProcessTask();

        treeWriter->Fill();

        modularDelphes->Clear();
        treeWriter->Clear();

        progressBar.Update(entry);
        ++entry;
      }
      progressBar.Finish();

      inputFile->Close();
    }

    modularDelphes->FinishTask();
    treeWriter->Write();

    cout << "** Exiting..." << endl;

    delete modularDelphes;
    delete confReader;
    delete treeWriter;
    delete outputFile;

    return 0;
  }
  catch(runtime_error &e)
  {
    if(treeWriter) delete treeWriter;
    if(outputFile) delete outputFile;
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}
