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

#include <iostream>
#include <fstream>
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
#include "modules/Delphes.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "HepMC3/ReaderAscii.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/GenPdfInfo.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Units.h"

using namespace std;
using namespace HepMC3;

map<Int_t, pair<Int_t, Int_t> > gMotherMap;
map<Int_t, pair<Int_t, Int_t> > gDaughterMap;

//---------------------------------------------------------------------------

void AnalyzeParticle(Bool_t in, Int_t counter,
  Double_t momentumCoefficient,
  Double_t positionCoefficient,
  shared_ptr<HepMC3::GenVertex> vertex,
  shared_ptr<HepMC3::GenParticle> particle,
  DelphesFactory *factory,
  TObjArray *allParticleOutputArray,
  TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray)
{
  map<Int_t, pair<Int_t, Int_t> >::iterator itMotherMap;
  map<Int_t, pair<Int_t, Int_t> >::iterator itDaughterMap;

  Candidate *candidate;
  TDatabasePDG *pdg;
  TParticlePDG *pdgParticle;
  Int_t pdgCode;

  Int_t pid, status, inVertexCode, outVertexCode;
  Double_t px, py, pz, e, mass;
  Double_t x, y, z, t;

  pdg = TDatabasePDG::Instance();

  candidate = factory->NewCandidate();

  pid = particle->pid();
  px = particle->momentum().px();
  py = particle->momentum().py();
  pz = particle->momentum().pz();
  e = particle->momentum().e();
  mass = particle->generated_mass();
  x = vertex->position().x();
  y = vertex->position().y();
  z = vertex->position().z();
  t = vertex->position().t();
  status = particle->status();

  outVertexCode = vertex->id();
  inVertexCode = particle->end_vertex() ? particle->end_vertex()->id() : 0;

  candidate->PID = pid;
  pdgCode = TMath::Abs(pid);

  candidate->Status = status;

  pdgParticle = pdg->GetParticle(pid);
  candidate->Charge = pdgParticle ? int(pdgParticle->Charge() / 3.0) : -999;
  candidate->Mass = mass;

  candidate->Momentum.SetPxPyPzE(px, py, pz, e);
  if(momentumCoefficient != 1.0)
  {
    candidate->Momentum *= momentumCoefficient;
  }

  candidate->M2 = 1;
  candidate->D2 = 1;

  if(in)
  {
    candidate->M1 = 1;
    candidate->Position.SetXYZT(0.0, 0.0, 0.0, 0.0);
  }
  else
  {
    candidate->M1 = outVertexCode;
    candidate->Position.SetXYZT(x, y, z, t);
    if(positionCoefficient != 1.0)
    {
      candidate->Position *= positionCoefficient;
    }

    itDaughterMap = gDaughterMap.find(outVertexCode);
    if(itDaughterMap == gDaughterMap.end())
    {
      gDaughterMap[outVertexCode] = make_pair(counter, counter);
    }
    else
    {
      itDaughterMap->second.second = counter;
    }
  }

  if(inVertexCode < 0)
  {
    candidate->D1 = inVertexCode;

    itMotherMap = gMotherMap.find(inVertexCode);
    if(itMotherMap == gMotherMap.end())
    {
      gMotherMap[inVertexCode] = make_pair(counter, -1);
    }
    else
    {
      itMotherMap->second.second = counter;
    }
  }
  else
  {
    candidate->D1 = 1;
  }

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

//---------------------------------------------------------------------------

void AnalyzeEvent(GenEvent &event, ExRootTreeBranch *branchEvent, DelphesFactory *factory,
  TObjArray *allParticleOutputArray, TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray, TStopwatch *readStopWatch, TStopwatch *procStopWatch)
{
  Int_t i, counter;
  map<Int_t, pair<Int_t, Int_t> >::iterator itMotherMap;
  map<Int_t, pair<Int_t, Int_t> >::iterator itDaughterMap;

  HepMCEvent *element;
  Candidate *candidate;
  Double_t momentumCoefficient, positionCoefficient;

  shared_ptr<IntAttribute> processID = event.attribute<IntAttribute>("signal_process_id");
  shared_ptr<IntAttribute> mpi = event.attribute<IntAttribute>("mpi");
  shared_ptr<DoubleAttribute> scale = event.attribute<DoubleAttribute>("event_scale");
  shared_ptr<DoubleAttribute> alphaQED = event.attribute<DoubleAttribute>("alphaQED");
  shared_ptr<DoubleAttribute> alphaQCD = event.attribute<DoubleAttribute>("alphaQCD");

  shared_ptr<GenCrossSection> cs = event.attribute<GenCrossSection>("GenCrossSection");
  shared_ptr<GenPdfInfo> pdf = event.attribute<GenPdfInfo>("GenPdfInfo");

  element = static_cast<HepMCEvent *>(branchEvent->NewEntry());

  element->Number = event.event_number();

  element->ProcessID = processID ? processID->value() : 0;
  element->MPI = mpi ? mpi->value() : 0;
  element->Weight = event.weights().size() > 0 ? event.weights()[0] : 1.0;
  element->Scale = scale ? scale->value() : 0.0;
  element->AlphaQED = alphaQED ? alphaQED->value() : 0.0;
  element->AlphaQCD = alphaQCD ? alphaQCD->value() : 0.0;

  if(cs)
  {
    element->CrossSection = cs->xsec();
    element->CrossSectionError = cs->xsec_err();
  }
  else
  {
    element->CrossSection = 0.0;
    element->CrossSectionError = 0.0;;
  }

  if(pdf)
  {
    element->ID1 = pdf->parton_id[0];
    element->ID2 = pdf->parton_id[1];
    element->X1 = pdf->x[0];
    element->X2 = pdf->x[1];
    element->ScalePDF = pdf->scale;
    element->PDF1 = pdf->xf[0];
    element->PDF2 = pdf->xf[1];
  }
  else
  {
    element->ID1 = 0;
    element->ID2 = 0;
    element->X1 = 0.0;
    element->X2 = 0.0;
    element->ScalePDF = 0.0;
    element->PDF1 = 0.0;
    element->PDF2 = 0.0;
  }

  element->ReadTime = readStopWatch->RealTime();
  element->ProcTime = procStopWatch->RealTime();

  momentumCoefficient = (event.momentum_unit() == Units::MEV) ? 0.001 : 1.0;
  positionCoefficient = (event.length_unit() == Units::MM) ? 1.0 : 10.0;

  counter = 0;
  for(auto vertex: event.vertices())
  {
    for(auto particle: vertex->particles_in())
    {
      if(!particle->production_vertex() || particle->production_vertex()->id() == 0)
      {
        AnalyzeParticle(kTRUE, counter, momentumCoefficient, positionCoefficient, vertex, particle,
          factory, allParticleOutputArray, stableParticleOutputArray, partonOutputArray);
        ++counter;
      }
    }
    for(auto particle: vertex->particles_out())
    {
      AnalyzeParticle(kFALSE, counter, momentumCoefficient, positionCoefficient, vertex, particle,
        factory, allParticleOutputArray, stableParticleOutputArray, partonOutputArray);
      ++counter;
    }
  }

  for(i = 0; i < allParticleOutputArray->GetEntriesFast(); ++i)
  {
    candidate = static_cast<Candidate *>(allParticleOutputArray->At(i));

    if(candidate->M1 > 0)
    {
      candidate->M1 = -1;
      candidate->M2 = -1;
    }
    else
    {
      itMotherMap = gMotherMap.find(candidate->M1);
      if(itMotherMap == gMotherMap.end())
      {
        candidate->M1 = -1;
        candidate->M2 = -1;
      }
      else
      {
        candidate->M1 = itMotherMap->second.first;
        candidate->M2 = itMotherMap->second.second;
      }
    }
    if(candidate->D1 > 0)
    {
      candidate->D1 = -1;
      candidate->D2 = -1;
    }
    else
    {
      itDaughterMap = gDaughterMap.find(candidate->D1);
      if(itDaughterMap == gDaughterMap.end())
      {
        candidate->D1 = -1;
        candidate->D2 = -1;
      }
      else
      {
        candidate->D1 = itDaughterMap->second.first;
        candidate->D2 = itDaughterMap->second.second;
      }
    }
  }
}

//---------------------------------------------------------------------------

void AnalyzeWeight(GenEvent &event, ExRootTreeBranch *branchWeight)
{
  Weight *element;

  for(auto weight: event.weights())
  {
    element = static_cast<Weight *>(branchWeight->NewEntry());

    element->Weight = weight;
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
  char appName[] = "DelphesHepMC3";
  stringstream message;
  ifstream inputFile;
  filebuf *inputBuffer;
  TFile *outputFile = 0;
  TStopwatch readStopWatch, procStopWatch;
  ExRootTreeWriter *treeWriter = 0;
  ExRootTreeBranch *branchEvent = 0, *branchWeight = 0;
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  DelphesFactory *factory = 0;
  TObjArray *allParticleOutputArray = 0, *stableParticleOutputArray = 0, *partonOutputArray = 0;
  ReaderAscii *reader = 0;
  GenEvent event;
  Int_t i, maxEvents, skipEvents;
  Long64_t length, eventCounter;

  if(argc < 3)
  {
    cout << " Usage: " << appName << " config_file"
         << " output_file"
         << " [input_file(s)]" << endl;
    cout << " config_file - configuration file in Tcl format," << endl;
    cout << " output_file - output file in ROOT format," << endl;
    cout << " input_file(s) - input file(s) in HepMC format," << endl;
    cout << " with no input_file, or when input_file is -, read standard input." << endl;
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
      message << "can't create output file " << argv[2];
      throw runtime_error(message.str());
    }

    treeWriter = new ExRootTreeWriter(outputFile, "Delphes");

    branchEvent = treeWriter->NewBranch("Event", HepMCEvent::Class());
    branchWeight = treeWriter->NewBranch("Weight", Weight::Class());

    confReader = new ExRootConfReader;
    confReader->ReadFile(argv[1]);

    maxEvents = confReader->GetInt("::MaxEvents", 0);
    skipEvents = confReader->GetInt("::SkipEvents", 0);

    if(maxEvents < 0)
    {
      throw runtime_error("MaxEvents must be zero or positive");
    }

    if(skipEvents < 0)
    {
      throw runtime_error("SkipEvents must be zero or positive");
    }

    modularDelphes = new Delphes("Delphes");
    modularDelphes->SetConfReader(confReader);
    modularDelphes->SetTreeWriter(treeWriter);

    factory = modularDelphes->GetFactory();
    allParticleOutputArray = modularDelphes->ExportArray("allParticles");
    stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
    partonOutputArray = modularDelphes->ExportArray("partons");

    inputBuffer = inputFile.rdbuf();
    reader = new ReaderAscii(inputFile);

    modularDelphes->InitTask();

    i = 3;
    do
    {
      if(interrupted) break;

      if(i == argc || strncmp(argv[i], "-", 2) == 0)
      {
        cout << "** Reading standard input" << endl;
        inputFile.ios::rdbuf(cin.rdbuf());
        length = -1;
      }
      else
      {
        cout << "** Reading " << argv[i] << endl;
        inputFile.ios::rdbuf(inputBuffer);
        inputFile.open(argv[i]);

        if(inputFile.fail())
        {
          message << "can't open " << argv[i];
          throw runtime_error(message.str());
        }

        inputFile.seekg(0, ios::end);
        length = inputFile.tellg();
        inputFile.seekg(0, ios::beg);

        if(length <= 0)
        {
          inputFile.close();
          inputFile.clear();
          ++i;
          continue;
        }
      }

      ExRootProgressBar progressBar(length);

      reader->skip(skipEvents);
      progressBar.Update(inputFile.tellg(), skipEvents);

      eventCounter = skipEvents;
      modularDelphes->Clear();
      treeWriter->Clear();
      event.clear();
      gMotherMap.clear();
      gDaughterMap.clear();
      readStopWatch.Start();
      reader->read_event(event);
      while((maxEvents <= 0 || eventCounter - skipEvents < maxEvents) && !reader->failed() && !interrupted)
      {
        readStopWatch.Stop();

        ++eventCounter;

        procStopWatch.Start();
        AnalyzeEvent(event, branchEvent, factory,
          allParticleOutputArray, stableParticleOutputArray,
          partonOutputArray, &readStopWatch, &procStopWatch);
        modularDelphes->ProcessTask();
        AnalyzeWeight(event, branchWeight);
        procStopWatch.Stop();

        treeWriter->Fill();

        modularDelphes->Clear();
        treeWriter->Clear();
        event.clear();
        gMotherMap.clear();
        gDaughterMap.clear();

        progressBar.Update(inputFile.tellg(), eventCounter);

        readStopWatch.Start();
        reader->read_event(event);
      }

      progressBar.Update(length, eventCounter, kTRUE);
      progressBar.Finish();

      if(length > 0) inputFile.close();

      ++i;
    } while(i < argc);

    modularDelphes->FinishTask();
    treeWriter->Write();

    cout << "** Exiting..." << endl;

    delete reader;
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
