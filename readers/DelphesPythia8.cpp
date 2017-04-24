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
#include <string>

#include <signal.h>

#include "Pythia.h"
#include "Pythia8Plugins/CombineMatchingInput.h"

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TObjArray.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

#include "modules/Delphes.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesLHEFReader.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

using namespace std;

//---------------------------------------------------------------------------

void ConvertInput(Long64_t eventCounter, Pythia8::Pythia *pythia,
  ExRootTreeBranch *branch, DelphesFactory *factory,
  TObjArray *allParticleOutputArray, TObjArray *stableParticleOutputArray, TObjArray *partonOutputArray,
  TStopwatch *readStopWatch, TStopwatch *procStopWatch)
{
  int i;

  HepMCEvent *element;
  Candidate *candidate;
  TDatabasePDG *pdg;
  TParticlePDG *pdgParticle;
  Int_t pdgCode;

  Int_t pid, status;
  Double_t px, py, pz, e, mass;
  Double_t x, y, z, t;

  // event information
  element = static_cast<HepMCEvent *>(branch->NewEntry());

  element->Number = eventCounter;

  element->ProcessID = pythia->info.code();
  element->MPI = 1;
  element->Weight = pythia->info.weight();
  element->Scale = pythia->info.QRen();
  element->AlphaQED = pythia->info.alphaEM();
  element->AlphaQCD = pythia->info.alphaS();

  element->ID1 = pythia->info.id1();
  element->ID2 = pythia->info.id2();
  element->X1 = pythia->info.x1();
  element->X2 = pythia->info.x2();
  element->ScalePDF = pythia->info.QFac();
  element->PDF1 = pythia->info.pdf1();
  element->PDF2 = pythia->info.pdf2();

  element->ReadTime = readStopWatch->RealTime();
  element->ProcTime = procStopWatch->RealTime();

  pdg = TDatabasePDG::Instance();

  for(i = 1; i < pythia->event.size(); ++i)
  {
    Pythia8::Particle &particle = pythia->event[i];

    pid = particle.id();
    status = particle.statusHepMC();
    px = particle.px(); py = particle.py(); pz = particle.pz(); e = particle.e(); mass = particle.m();
    x = particle.xProd(); y = particle.yProd(); z = particle.zProd(); t = particle.tProd();

    candidate = factory->NewCandidate();

    candidate->PID = pid;
    pdgCode = TMath::Abs(candidate->PID);

    candidate->Status = status;

    candidate->M1 = particle.mother1() - 1;
    candidate->M2 = particle.mother2() - 1;

    candidate->D1 = particle.daughter1() - 1;
    candidate->D2 = particle.daughter2() - 1;

    pdgParticle = pdg->GetParticle(pid);
    candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
    candidate->Mass = mass;

    candidate->Momentum.SetPxPyPzE(px, py, pz, e);

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

/*
Single-particle gun. The particle must be a colour singlet.
Input: flavour, energy, direction (theta, phi).
If theta < 0 then random choice over solid angle.
Optional final argument to put particle at rest => E = m.
from pythia8 example 21
*/

void fillParticle(int id, double pMax, double etaMax,
  Pythia8::Event &event, Pythia8::ParticleData &pdt, Pythia8::Rndm &rndm)
{
  // Reset event record to allow for new event.
  event.reset();

  // Generate uniform pt and eta.
  double pt, eta, phi, pp, ee, mm;

  // pMin = 0.1 GeV for single particles
  pp = pow(10, - 1.0 + (log10(pMax) + 1.0) * rndm.flat());
  eta = (2.0 * rndm.flat() - 1.0) * etaMax;
  phi = 2.0 * M_PI * rndm.flat();
  mm = pdt.mSel(id);
  ee = Pythia8::sqrtpos(pp*pp + mm*mm);
  pt = pp / cosh(eta);

  // Store the particle in the event record.
  event.append(id, 1, 0, 0, pt * cos(phi), pt * sin(phi), pt * sinh(eta), ee, mm);
}

//---------------------------------------------------------------------------

void fillPartons(int id, double pMax, double etaMax,
  Pythia8::Event &event, Pythia8::ParticleData &pdt, Pythia8::Rndm &rndm)
{
  // Reset event record to allow for new event.
  event.reset();

  // Generate uniform pt and eta.
  double pt, eta, phi, pp, ee, mm;

  // pMin = 1 GeV for jets
  pp = pow(10, log10(pMax) * rndm.flat());
  eta = (2.0 * rndm.flat() - 1.0) * etaMax;
  phi = 2.0 * M_PI * rndm.flat();
  mm = pdt.mSel(id);
  ee = Pythia8::sqrtpos(pp*pp + mm*mm);
  pt = pp / cosh(eta);

  if( (id == 4 || id == 5) && pt < 10.0) return;

  if(id == 21)
  {
    event.append(21, 23, 101, 102, pt * cos(phi), pt * sin(phi), pt * sinh(eta), ee);
    event.append(21, 23, 102, 101, -pt * cos(phi), -pt * sin(phi), -pt * sinh(eta), ee);
  }
  else
  {
    event.append(id, 23, 101, 0, pt * cos(phi), pt * sin(phi), pt * sinh(eta), ee, mm);
    event.append(-id, 23, 0, 101, -pt * cos(phi), -pt * sin(phi), -pt * sinh(eta), ee, mm);
  }
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "DelphesPythia8";
  stringstream message;
  FILE *inputFile = 0;
  TFile *outputFile = 0;
  TStopwatch readStopWatch, procStopWatch;
  ExRootTreeWriter *treeWriter = 0;
  ExRootTreeBranch *branchEvent = 0;
  ExRootTreeBranch *branchEventLHEF = 0, *branchWeightLHEF = 0;
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  DelphesFactory *factory = 0;
  TObjArray *stableParticleOutputArray = 0, *allParticleOutputArray = 0, *partonOutputArray = 0;
  TObjArray *stableParticleOutputArrayLHEF = 0, *allParticleOutputArrayLHEF = 0, *partonOutputArrayLHEF = 0;
  DelphesLHEFReader *reader = 0;
  Long64_t eventCounter, errorCounter;
  Long64_t numberOfEvents, timesAllowErrors;
  Bool_t spareFlag1;
  Int_t spareMode1;
  Double_t spareParm1, spareParm2;

  Pythia8::Pythia *pythia = 0;

  // for matching
  Pythia8::CombineMatchingInput *combined = 0;
  Pythia8::UserHooks* matching = 0;

  if(argc != 4)
  {
    cout << " Usage: " << appName << " config_file" << " pythia_card" << " output_file" << endl;
    cout << " config_file - configuration file in Tcl format," << endl;
    cout << " pythia_card - Pythia8 configuration file," << endl;
    cout << " output_file - output file in ROOT format." << endl;
    return 1;
  }

  signal(SIGINT, SignalHandler);

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  try
  {
    outputFile = TFile::Open(argv[3], "CREATE");

    if(outputFile == NULL)
    {
      message << "can't create output file " << argv[3];
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

    // Initialize Pythia
    pythia = new Pythia8::Pythia;

    // jet matching
    matching = combined->getHook(*pythia);
    if(!matching)
    {
      throw runtime_error("can't do matching");
    }
    pythia->setUserHooksPtr(matching);


    if(pythia == NULL)
    {
      throw runtime_error("can't create Pythia instance");
    }

    // Read in commands from configuration file
    if(!pythia->readFile(argv[2]))
    {
      message << "can't read Pythia8 configuration file " << argv[2] << endl;
      throw runtime_error(message.str());
    }

    // Extract settings to be used in the main program
    numberOfEvents = pythia->mode("Main:numberOfEvents");
    timesAllowErrors = pythia->mode("Main:timesAllowErrors");

    spareFlag1 = pythia->flag("Main:spareFlag1");
    spareMode1 = pythia->mode("Main:spareMode1");
    spareParm1 = pythia->parm("Main:spareParm1");
    spareParm2 = pythia->parm("Main:spareParm2");

    // Check if particle gun
    if(!spareFlag1)
    {
      inputFile = fopen(pythia->word("Beams:LHEF").c_str(), "r");
      if(inputFile)
      {
        reader = new DelphesLHEFReader;
        reader->SetInputFile(inputFile);

        branchEventLHEF = treeWriter->NewBranch("EventLHEF", LHEFEvent::Class());
        branchWeightLHEF = treeWriter->NewBranch("WeightLHEF", LHEFWeight::Class());

        allParticleOutputArrayLHEF = modularDelphes->ExportArray("allParticlesLHEF");
        stableParticleOutputArrayLHEF = modularDelphes->ExportArray("stableParticlesLHEF");
        partonOutputArrayLHEF = modularDelphes->ExportArray("partonsLHEF");
      }
    }

    modularDelphes->InitTask();

    pythia->init();

    // ExRootProgressBar progressBar(numberOfEvents - 1);
    ExRootProgressBar progressBar(-1);

    // Loop over all events
    errorCounter = 0;
    treeWriter->Clear();
    modularDelphes->Clear();
    readStopWatch.Start();
    for(eventCounter = 0; eventCounter < numberOfEvents && !interrupted; ++eventCounter)
    {
      while(reader && reader->ReadBlock(factory, allParticleOutputArrayLHEF,
        stableParticleOutputArrayLHEF, partonOutputArrayLHEF) && !reader->EventReady());

      if(spareFlag1)
      {
        if((spareMode1 >= 1 && spareMode1 <= 5) || spareMode1 == 21)
        {
          fillPartons(spareMode1, spareParm1, spareParm2, pythia->event, pythia->particleData, pythia->rndm);
        }
        else
        {
          fillParticle(spareMode1, spareParm1, spareParm2, pythia->event, pythia->particleData, pythia->rndm);
        }
      }

      if(!pythia->next())
      {
        // If failure because reached end of file then exit event loop
        if(pythia->info.atEndOfFile())
        {
          cerr << "Aborted since reached end of Les Houches Event File" << endl;
          break;
        }

        // First few failures write off as "acceptable" errors, then quit
        if(++errorCounter > timesAllowErrors)
        {
          cerr << "Event generation aborted prematurely, owing to error!" << endl;
          break;
        }

        modularDelphes->Clear();
        reader->Clear();
        continue;
      }

      readStopWatch.Stop();

      procStopWatch.Start();
      ConvertInput(eventCounter, pythia, branchEvent, factory,
        allParticleOutputArray, stableParticleOutputArray, partonOutputArray,
        &readStopWatch, &procStopWatch);
      modularDelphes->ProcessTask();
      procStopWatch.Stop();

      if(reader)
      {
        reader->AnalyzeEvent(branchEventLHEF, eventCounter, &readStopWatch, &procStopWatch);
        reader->AnalyzeWeight(branchWeightLHEF);
      }

      treeWriter->Fill();

      treeWriter->Clear();
      modularDelphes->Clear();
      if(reader) reader->Clear();

      readStopWatch.Start();
      progressBar.Update(eventCounter, eventCounter);
    }

    progressBar.Update(eventCounter, eventCounter, kTRUE);
    progressBar.Finish();

    pythia->stat();

    modularDelphes->FinishTask();
    treeWriter->Write();

    cout << "** Exiting..." << endl;

    delete reader;
    delete pythia;
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
