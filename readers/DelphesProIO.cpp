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
 *  @author S.Chekanov (ANL)
 */

#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include <map>

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

#include <proio/event.h>
#include <proio/model/mc.pb.h>
#include <proio/reader.h>
namespace model = proio::model::mc;

using namespace std;

//---------------------------------------------------------------------------
// This method dynamically checks the message type (varint or not) depending on
// non-zero value of units momentumUnit and positionUnit.
void ConvertInput(proio::Event *event, double momentumUnit, double positionUnit,
  ExRootTreeBranch *branch, DelphesFactory *factory,
  TObjArray *allParticleOutputArray, TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray, TStopwatch *readStopWatch, TStopwatch *procStopWatch)
{

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
  element = static_cast<HepMCEvent *>(branch->NewEntry());

  int nID = 0;
  double weight = 0;
  int process_id = 0;
  auto mcentries = event->TaggedEntries("MCParameters");
  for(uint64_t mcentryID : mcentries)
  {
    auto mcpar = dynamic_cast<proio::model::mc::MCParameters *>(event->GetEntry(mcentryID));
    nID = mcpar->number();
    weight = mcpar->weight();
    process_id = mcpar->processid();
    break; // consider only most generic from 1st entry
  };

  element->Number = nID;
  element->ProcessID = process_id;
  element->Weight = weight;

  /*
  // Pythia8 specific
  element->MPI = mutableEvent->mpi();
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
*/

  element->ReadTime = readStopWatch->RealTime();
  element->ProcTime = procStopWatch->RealTime();

  if(momentumUnit > 0 && positionUnit > 0)
  {

    auto entries = event->TaggedEntries("VarintPackedParticles");

    for(uint64_t entryID : entries)
    {

      auto mutableParticles = dynamic_cast<model::VarintPackedParticles *>(event->GetEntry(entryID));

      for(int i = 0; i < mutableParticles->pdg_size(); ++i)
      {
        pid = mutableParticles->pdg(i);
        status = mutableParticles->status(i);
        px = mutableParticles->px(i) / momentumUnit;
        py = mutableParticles->py(i) / momentumUnit;
        pz = mutableParticles->pz(i) / momentumUnit;
        mass = mutableParticles->mass(i) / momentumUnit;
        x = mutableParticles->x(i) / positionUnit;
        y = mutableParticles->y(i) / positionUnit;
        z = mutableParticles->z(i) / positionUnit;
        t = mutableParticles->t(i) / positionUnit;

        candidate = factory->NewCandidate();
        candidate->PID = pid;
        pdgCode = TMath::Abs(candidate->PID);
        candidate->Status = status;
        candidate->M1 = mutableParticles->parent1(i);
        candidate->M2 = mutableParticles->parent2(i);
        candidate->D1 = mutableParticles->child1(i);
        candidate->D2 = mutableParticles->child2(i);
        pdgParticle = pdg->GetParticle(pid);
        candidate->Charge = mutableParticles->charge(i) / 3.0;
        //candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
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
  }
  else
  {

    auto entries = event->TaggedEntries("Particle");

    for(uint64_t entryID : entries)
    {
      auto part = dynamic_cast<model::Particle *>(event->GetEntry(entryID));
      pid = part->pdg();
      status = part->status();
      model::XYZF pvector = part->p();
      px = pvector.x();
      py = pvector.y();
      pz = pvector.z();
      mass = part->mass();
      auto v = part->vertex();
      x = v.x();
      y = v.y();
      z = v.z();
      t = v.t();

      candidate = factory->NewCandidate();
      candidate->PID = pid;
      pdgCode = TMath::Abs(candidate->PID);
      candidate->Status = status;

      int M1 = 0;
      int M2 = 0;
      for(int k1 = 0; k1 < part->parent_size(); k1++)
      {
        if(k1 == 0)
        {
          auto mother = dynamic_cast<model::Particle *>(event->GetEntry(part->parent(0)));
          M1 = mother->barcode();
        }
        if(k1 == 1)
        {
          auto mother = dynamic_cast<model::Particle *>(event->GetEntry(part->parent(1)));
          M2 = mother->barcode();
        }
      }

      int D1 = 0;
      int D2 = 0;
      for(int k1 = 0; k1 < part->child_size(); k1++)
      {
        if(k1 == 0)
        {
          auto child = dynamic_cast<model::Particle *>(event->GetEntry(part->child(0)));
          D1 = child->barcode();
        }

        if(k1 == 1)
        {
          auto child = dynamic_cast<model::Particle *>(event->GetEntry(part->child(1)));
          D2 = child->barcode();
        };
      }

      candidate->M1 = M1;
      candidate->M2 = M2;
      candidate->D1 = D1;
      candidate->D2 = D2;

      pdgParticle = pdg->GetParticle(pid);
      candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge() / 3.0) : -999;
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

  } // end particle type
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
  char appName[] = "DelphesProIO";
  stringstream message;
  proio::Reader *inputFile = 0;
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

  if(argc < 4)
  {
    cout << " Usage: " << appName << " config_file"
         << " output_file"
         << " input_file(s)" << endl;
    cout << " config_file - configuration file in Tcl format," << endl;
    cout << " output_file - output file in ROOT format," << endl;
    cout << " input_file(s) - input file(s) in ProIO format." << endl;
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
      cout << "** INFO: Reading " << argv[i] << endl;

      inputFile = new proio::Reader(argv[i]);

      if(inputFile == NULL)
      {
        message << "can't open " << argv[i] << endl;
        throw runtime_error(message.str());
      }

      /*
// this is slow method, but general
      inputFile->SeekToStart(); 
      int nn=0;
      auto event = new proio::Event();
      while(inputFile->Next(event)){
            auto entries = event->TaggedEntries("Particle");
            if (entries.size()>0) nn++;
      }
*/

      auto event = new proio::Event();

      double varint_energy = 0;
      double varint_length = 0;

      auto max_n_events = std::numeric_limits<uint64_t>::max();
      auto nn = inputFile->Skip(max_n_events);
      cout << "** INFO: " << nn - 1 << " events found in ProIO file" << endl;
      inputFile->SeekToStart();
      numberOfEvents = nn - 1; // last event has only metadata (no particle record)

      if(numberOfEvents <= 0) continue;

      ExRootProgressBar progressBar(numberOfEvents - 1);

      // Loop over all objects
      modularDelphes->Clear();
      treeWriter->Clear();
      readStopWatch.Start();

      for(eventCounter = 0; eventCounter < numberOfEvents && !interrupted; ++eventCounter)
      {
        inputFile->Next(event);
        if(event == 0) continue;

        // get metadata
        if(eventCounter == 0)
        {
          auto metadata = event->Metadata();
          std::cout << "** INFO: ProIO file metadata:" << std::endl;
          for(auto element : metadata)
          {
            string key = (string)element.first;
            string value = (string)(*element.second);
            std::cout << "** INFO:   " << key << " = " << value << std::endl;
            if(key == "info:varint_energy") varint_energy = std::stod(value);
            if(key == "info:varint_length") varint_length = std::stod(value);
          }
        }

        readStopWatch.Stop();

        procStopWatch.Start();

        ConvertInput(event, varint_energy, varint_length,
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

      delete inputFile;
    }

    modularDelphes->FinishTask();
    treeWriter->Write();

    cout << "** Exiting..." << endl;

    delete modularDelphes;
    delete confReader;
    delete treeWriter;
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
