/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include <map>
#include <vector>

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

#include "TApplication.h"
#include "TROOT.h"

#include "TDatabasePDG.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TParticlePDG.h"
#include "TStopwatch.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesStream.h"
#include "modules/Delphes.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/WeightsInfo.h"

using namespace std;

//---------------------------------------------------------------------------

void ConvertInput(fwlite::Event &event, Long64_t eventCounter,
  ExRootTreeBranch *branchEvent, ExRootTreeBranch *branchWeight,
  DelphesFactory *factory, TObjArray *allParticleOutputArray,
  TObjArray *stableParticleOutputArray, TObjArray *partonOutputArray, Bool_t firstEvent)
{

  fwlite::Handle<GenEventInfoProduct> handleGenEventInfo;
  fwlite::Handle<LHEEventProduct> handleLHEEvent;
  fwlite::Handle<vector<reco::GenParticle> > handleParticle;
  fwlite::Handle<vector<pat::PackedGenParticle> > handlePackedParticle;

  vector<reco::GenParticle>::const_iterator itParticle;
  vector<pat::PackedGenParticle>::const_iterator itPackedParticle;

  vector<const reco::Candidate *> vectorCandidate;
  vector<const reco::Candidate *>::iterator itCandidate;

  handleGenEventInfo.getByLabel(event, "generator");

  if(!((handleLHEEvent.getBranchNameFor(event, "source")).empty()))
  {
    handleLHEEvent.getByLabel(event, "source");
  }
  else if(!((handleLHEEvent.getBranchNameFor(event, "externalLHEProducer")).empty()))
  {
    handleLHEEvent.getByLabel(event, "externalLHEProducer");
  }
  else if(firstEvent)
  {
    cout << "Wrong LHEEvent Label! Please, check the input file." << endl;
  }

  if(!((handleParticle.getBranchNameFor(event, "genParticles")).empty()))
  {
    handleParticle.getByLabel(event, "genParticles");
  }
  else if(!((handlePackedParticle.getBranchNameFor(event, "packedGenParticles")).empty()) && !((handleParticle.getBranchNameFor(event, "prunedGenParticles")).empty()))
  {
    handleParticle.getByLabel(event, "prunedGenParticles");
    handlePackedParticle.getByLabel(event, "packedGenParticles");
  }
  else
  {
    std::cout << "Wrong GenParticle Label! Please, check the input file." << std::endl;
    exit(-1);
  }

  Bool_t foundLHE = !((handleLHEEvent.getBranchNameFor(event, "source")).empty()) || !((handleLHEEvent.getBranchNameFor(event, "externalLHEProducer")).empty());
  Bool_t isMiniAOD = !((handlePackedParticle.getBranchNameFor(event, "packedGenParticles")).empty()) && ((handleParticle.getBranchNameFor(event, "genParticles")).empty());

  HepMCEvent *element;
  Weight *weight;
  Candidate *candidate;
  TDatabasePDG *pdg;
  TParticlePDG *pdgParticle;
  Int_t pdgCode;

  Int_t pid, status;
  Double_t px, py, pz, e, mass;
  Double_t x, y, z;

  element = static_cast<HepMCEvent *>(branchEvent->NewEntry());

  element->Number = eventCounter;

  element->ProcessID = handleGenEventInfo->signalProcessID();
  element->MPI = 1;
  element->Weight = handleGenEventInfo->weight();
  element->Scale = handleGenEventInfo->qScale();
  element->AlphaQED = handleGenEventInfo->alphaQED();
  element->AlphaQCD = handleGenEventInfo->alphaQCD();

  element->ID1 = 0;
  element->ID2 = 0;
  element->X1 = 0.0;
  element->X2 = 0.0;
  element->ScalePDF = 0.0;
  element->PDF1 = 0.0;
  element->PDF2 = 0.0;

  element->ReadTime = 0.0;
  element->ProcTime = 0.0;

  if(foundLHE)
  {
    const vector<gen::WeightsInfo> &vectorWeightsInfo = handleLHEEvent->weights();
    vector<gen::WeightsInfo>::const_iterator itWeightsInfo;

    for(itWeightsInfo = vectorWeightsInfo.begin(); itWeightsInfo != vectorWeightsInfo.end(); ++itWeightsInfo)
    {
      weight = static_cast<Weight *>(branchWeight->NewEntry());
      weight->Weight = itWeightsInfo->wgt;
    }
  }

  pdg = TDatabasePDG::Instance();

  for(itParticle = handleParticle->begin(); itParticle != handleParticle->end(); ++itParticle)
  {
    const reco::GenParticle &particle = *itParticle;
    if(!isMiniAOD || particle.status() != 1) vectorCandidate.push_back(&*itParticle);
  }

  for(itParticle = handleParticle->begin(); itParticle != handleParticle->end(); ++itParticle)
  {
    const reco::GenParticle &particle = *itParticle;

    pid = particle.pdgId();
    status = particle.status();
    if(isMiniAOD && particle.status() == 1) continue;
    px = particle.px();
    py = particle.py();
    pz = particle.pz();
    e = particle.energy();
    mass = particle.mass();
    x = particle.vx();
    y = particle.vy();
    z = particle.vz();

    candidate = factory->NewCandidate();

    candidate->PID = pid;
    pdgCode = TMath::Abs(candidate->PID);

    candidate->Status = status;

    if(particle.mother())
    {
      itCandidate = find(vectorCandidate.begin(), vectorCandidate.end(), particle.mother());
      if(itCandidate != vectorCandidate.end()) candidate->M1 = distance(vectorCandidate.begin(), itCandidate);
    }

    itCandidate = find(vectorCandidate.begin(), vectorCandidate.end(), particle.daughter(0));
    if(itCandidate != vectorCandidate.end()) candidate->D1 = distance(vectorCandidate.begin(), itCandidate);

    itCandidate = find(vectorCandidate.begin(), vectorCandidate.end(), particle.daughter(particle.numberOfDaughters() - 1));
    if(itCandidate != vectorCandidate.end()) candidate->D2 = distance(vectorCandidate.begin(), itCandidate);

    pdgParticle = pdg->GetParticle(pid);
    candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge() / 3.0) : -999;
    candidate->Mass = mass;

    candidate->Momentum.SetPxPyPzE(px, py, pz, e);

    candidate->Position.SetXYZT(x * 10.0, y * 10.0, z * 10.0, 0.0);

    allParticleOutputArray->Add(candidate);

    if(!pdgParticle) continue;

    if(status == 1)
    {
      // Prevent duplicated particle.
      if(!isMiniAOD) stableParticleOutputArray->Add(candidate);
      if (pdgCode == 11 || pdgCode == 13) partonOutputArray->Add(candidate);
    }
    //else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15)
    else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 11 || pdgCode == 13 || pdgCode == 15)
    {
      partonOutputArray->Add(candidate);
    }
  }

  if(!isMiniAOD) return;
  // For MiniAOD sample,
  // Only status==1 particles are saved to packedGenParticles.
  for(itPackedParticle = handlePackedParticle->begin(); itPackedParticle != handlePackedParticle->end(); ++itPackedParticle)
  {
    vectorCandidate.push_back(&*itPackedParticle);
  }

  for(itPackedParticle = handlePackedParticle->begin(); itPackedParticle != handlePackedParticle->end(); ++itPackedParticle)
  {
    const pat::PackedGenParticle &particle = *itPackedParticle;

    pid = particle.pdgId();
    status = particle.status();
    px = particle.px();
    py = particle.py();
    pz = particle.pz();
    e = particle.energy();
    mass = particle.mass();
    x = particle.vx();
    y = particle.vy();
    z = particle.vz();

    candidate = factory->NewCandidate();

    candidate->PID = pid;
    pdgCode = TMath::Abs(candidate->PID);

    candidate->Status = status;

    if(particle.mother(0))
    {
      itCandidate = find(vectorCandidate.begin(), vectorCandidate.end(), particle.mother(0));
      if(itCandidate != vectorCandidate.end()) candidate->M1 = distance(vectorCandidate.begin(), itCandidate);
    }

    itCandidate = find(vectorCandidate.begin(), vectorCandidate.end(), particle.daughter(0));
    if(itCandidate != vectorCandidate.end()) candidate->D1 = distance(vectorCandidate.begin(), itCandidate);

    itCandidate = find(vectorCandidate.begin(), vectorCandidate.end(), particle.daughter(particle.numberOfDaughters() - 1));
    if(itCandidate != vectorCandidate.end()) candidate->D2 = distance(vectorCandidate.begin(), itCandidate);

    pdgParticle = pdg->GetParticle(pid);
    candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge() / 3.0) : -999;
    candidate->Mass = mass;

    candidate->Momentum.SetPxPyPzE(px, py, pz, e);

    candidate->Position.SetXYZT(x * 10.0, y * 10.0, z * 10.0, 0.0);

    allParticleOutputArray->Add(candidate);

    if(!pdgParticle) continue;

    if(status == 1)
    {
      stableParticleOutputArray->Add(candidate);
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
  ExRootTreeBranch *branchEvent = 0, *branchWeight = 0;
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  DelphesFactory *factory = 0;
  TObjArray *allParticleOutputArray = 0, *stableParticleOutputArray = 0, *partonOutputArray = 0;
  Int_t i;
  Long64_t eventCounter, numberOfEvents;
  Bool_t firstEvent = kTRUE;

  if(argc < 4)
  {
    cout << " Usage: " << appName << " config_file"
         << " output_file"
         << " input_file(s)" << endl;
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

  FWLiteEnabler::enable();

  try
  {
    outputFile = TFile::Open(argv[2], "CREATE");

    if(outputFile == NULL)
    {
      message << "can't open " << argv[2] << endl;
      throw runtime_error(message.str());
    }

    treeWriter = new ExRootTreeWriter(outputFile, "Delphes");

    branchEvent = treeWriter->NewBranch("Event", HepMCEvent::Class());
    branchWeight = treeWriter->NewBranch("Weight", Weight::Class());

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

      numberOfEvents = event.size();

      if(numberOfEvents <= 0) continue;

      // ExRootProgressBar progressBar(numberOfEvents - 1);
      ExRootProgressBar progressBar(-1);

      // Loop over all objects
      eventCounter = 0;
      modularDelphes->Clear();
      treeWriter->Clear();

      for(event.toBegin(); !event.atEnd() && !interrupted; ++event)
      {
        ConvertInput(event, eventCounter, branchEvent, branchWeight, factory,
          allParticleOutputArray, stableParticleOutputArray, partonOutputArray, firstEvent);
        modularDelphes->ProcessTask();

        firstEvent = kFALSE;

        treeWriter->Fill();

        modularDelphes->Clear();
        treeWriter->Clear();

        progressBar.Update(eventCounter, eventCounter);
        ++eventCounter;
      }

      progressBar.Update(eventCounter, eventCounter, kTRUE);
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
