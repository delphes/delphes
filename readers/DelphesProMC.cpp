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

#include "ProMC.pb.h"
#include "ProMCBook.h"
#include "ProMCHeader.pb.h"

using namespace std;

//---------------------------------------------------------------------------

void ConvertInput(ProMCEvent &event, double momentumUnit, double positionUnit,
  ExRootTreeBranch *branch, DelphesFactory *factory,
  TObjArray *allParticleOutputArray, TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray, TStopwatch *readStopWatch, TStopwatch *procStopWatch)
{
  Int_t i;

  ProMCEvent_Event *mutableEvent;
  ProMCEvent_Particles *mutableParticles;

  HepMCEvent *element;
  Candidate *candidate;
  TDatabasePDG *pdg;
  TParticlePDG *pdgParticle;
  Int_t pdgCode;

  Int_t pid, status;
  Double_t px, py, pz, mass;
  Double_t x, y, z, t;

  pdg = TDatabasePDG::Instance();

  // event information
  mutableEvent = event.mutable_event();

  element = static_cast<HepMCEvent *>(branch->NewEntry());

  element->Number = mutableEvent->number();

  element->ProcessID = mutableEvent->process_id();
  element->MPI = mutableEvent->mpi();
  element->Weight = mutableEvent->weight();
  element->Scale = mutableEvent->scale();
  element->AlphaQED = mutableEvent->alpha_qed();
  element->AlphaQCD = mutableEvent->alpha_qcd();

  element->ID1 = mutableEvent->id1();
  element->ID2 = mutableEvent->id2();
  element->X1 = mutableEvent->x1();
  element->X2 = mutableEvent->x2();
  element->ScalePDF = mutableEvent->scale_pdf();
  element->PDF1 = mutableEvent->pdf1();
  element->PDF2 = mutableEvent->pdf2();

  element->ReadTime = readStopWatch->RealTime();
  element->ProcTime = procStopWatch->RealTime();

  mutableParticles = event.mutable_particles();

  for(i = 0; i < mutableParticles->pdg_id_size(); ++i)
  {
    pid = mutableParticles->pdg_id(i);
    status = mutableParticles->status(i);

    px = mutableParticles->px(i)/momentumUnit;
    py = mutableParticles->py(i)/momentumUnit;
    pz = mutableParticles->pz(i)/momentumUnit;
    mass = mutableParticles->mass(i)/momentumUnit;
    x = mutableParticles->x(i)/positionUnit;
    y = mutableParticles->y(i)/positionUnit;
    z = mutableParticles->z(i)/positionUnit;
    t = mutableParticles->t(i)/positionUnit;

    candidate = factory->NewCandidate();

    candidate->PID = pid;
    pdgCode = TMath::Abs(candidate->PID);

    candidate->Status = status;

    candidate->M1 = mutableParticles->mother1(i);
    candidate->M2 = mutableParticles->mother2(i);

    candidate->D1 = mutableParticles->daughter1(i);
    candidate->D2 = mutableParticles->daughter2(i);

    pdgParticle = pdg->GetParticle(pid);
    candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
    candidate->Mass = mass;

    candidate->Momentum.SetXYZM(px, py, pz, mass);

    candidate->Position.SetXYZT(x, y, z, t);

    allParticleOutputArray->Add(candidate);

    if(!pdgParticle) continue;

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

//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "DelphesProMC";
  stringstream message;
  ProMCBook *inputFile = 0;
  TFile *outputFile = 0;
  TStopwatch readStopWatch, procStopWatch;
  ExRootTreeWriter *treeWriter = 0;
  ExRootTreeBranch *branchEvent = 0;
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  DelphesFactory *factory = 0;
  TObjArray *allParticleOutputArray = 0, *stableParticleOutputArray = 0, *partonOutputArray = 0;
  Int_t i;
  Long64_t eventCounter, numberOfEvents;
  double momentumUnit = 1.0, positionUnit = 1.0;

  if(argc < 4)
  {
    cout << " Usage: " << appName << " config_file" << " output_file" << " input_file(s)" << endl;
    cout << " config_file - configuration file in Tcl format," << endl;
    cout << " output_file - output file in ROOT format," << endl;
    cout << " input_file(s) - input file(s) in ProMC format." << endl;
    return 1;
  }

  signal(SIGINT, SignalHandler);

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

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

      inputFile = new ProMCBook(argv[i], "r");

      ProMCHeader header = inputFile->getHeader();

      momentumUnit = static_cast<double>(header.momentumunit());
      positionUnit = static_cast<double>(header.lengthunit());



      if(inputFile == NULL)
      {
        message << "can't open " << argv[i] << endl;
        throw runtime_error(message.str());
      }

      numberOfEvents = inputFile->getEvents();

      if(numberOfEvents <= 0) continue;

      ExRootProgressBar progressBar(numberOfEvents - 1);

      // Loop over all objects
      modularDelphes->Clear();
      treeWriter->Clear();
      readStopWatch.Start();
      for(eventCounter = 0; eventCounter < numberOfEvents && !interrupted; ++eventCounter)
      {
        if(inputFile->next() != 0) continue;
        ProMCEvent event = inputFile->get();

        readStopWatch.Stop();

        procStopWatch.Start();
        ConvertInput(event, momentumUnit, positionUnit,
          branchEvent, factory,
          allParticleOutputArray, stableParticleOutputArray,
          partonOutputArray, &readStopWatch, &procStopWatch);
        modularDelphes->ProcessTask();
        procStopWatch.Stop();

        treeWriter->Fill();

        modularDelphes->Clear();
        treeWriter->Clear();

        readStopWatch.Start();
        progressBar.Update(eventCounter);
      }

      progressBar.Update(eventCounter, eventCounter, kTRUE);
      progressBar.Finish();

      inputFile->close();
      delete inputFile;
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
