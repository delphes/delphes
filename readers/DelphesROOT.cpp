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
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <memory>

#include <map>
#include <vector>

#include <stdlib.h>
#include <signal.h>
#include <stdio.h>

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TClonesArray.h"
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
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"



using namespace std;

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "DelphesROOT";
  stringstream message;
  TFile *inputFile = 0;
  TFile *outputFile = 0;
  TStopwatch eventStopWatch;
  ExRootTreeWriter *treeWriter = 0;
  ExRootTreeBranch *branchEvent = 0;
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  DelphesFactory *factory = 0;
  GenParticle *gen;
  HepMCEvent *element, *eve;
  Candidate *candidate;
  Int_t pdgCode;
  Int_t maxEvents;

  const Double_t c_light = 2.99792458E8;

  TObjArray *allParticleOutputArray = 0, *stableParticleOutputArray = 0, *partonOutputArray = 0;
  Int_t i;
  Long64_t eventCounter, numberOfEvents;

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
   
    maxEvents = confReader->GetInt("::MaxEvents", 0);

    if(maxEvents < 0)
    {
      throw runtime_error("MaxEvents must be zero or positive");
    }

    modularDelphes = new Delphes("Delphes");
    modularDelphes->SetConfReader(confReader);
    modularDelphes->SetTreeWriter(treeWriter);
    
    TChain *chain = new TChain("Delphes");
    
    factory = modularDelphes->GetFactory();
    allParticleOutputArray = modularDelphes->ExportArray("allParticles");
    stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
    partonOutputArray = modularDelphes->ExportArray("partons");

    modularDelphes->InitTask();

    for(i = 3; i < argc && !interrupted; ++i)
    {
      cout << "** Reading " << argv[i] << endl;

      chain->Add(argv[i]);
      ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
      
      inputFile = TFile::Open(argv[i]);

      if(inputFile == NULL)
      {
        message << "can't open " << argv[i] << endl;
        throw runtime_error(message.str());
      }
      
      numberOfEvents = treeReader->GetEntries();
      TClonesArray *branchParticle   = treeReader->UseBranch("Particle");
      TClonesArray *branchHepMCEvent = treeReader->UseBranch("Event");
     
      if(numberOfEvents <= 0) continue;

      // ExRootProgressBar progressBar(numberOfEvents - 1);
      ExRootProgressBar progressBar(-1);

      // Loop over all objects
      eventCounter = 0;
      modularDelphes->Clear();
      treeWriter->Clear();
      for(Int_t entry = 0; entry < numberOfEvents && !interrupted; ++entry)
      {
    
        if(entry >= maxEvents) break;
        treeReader->ReadEntry(entry);
       
        // -- TBC need also to include event weights --  
        
        eve = (HepMCEvent*) branchHepMCEvent->At(0);
        element = static_cast<HepMCEvent *>(branchEvent->NewEntry());
        
        element->Number = eventCounter;

        element->ProcessID = eve->ProcessID;
        element->MPI = eve->MPI;
        element->Weight = eve->Weight;
        element->Scale = eve->Scale;
        element->AlphaQED = eve->AlphaQED;
        element->AlphaQCD = eve->AlphaQCD;

        element->ID1 = eve->ID1;
        element->ID2 = eve->ID2;
        element->X1 = eve->X1;
        element->X2 = eve->X2;
        element->ScalePDF = eve->ScalePDF;
        element->PDF1 = eve->PDF1;
        element->PDF2 = eve->PDF2;

        element->ReadTime = eve->ReadTime;
        element->ProcTime = eve->ProcTime;

        for(Int_t j=0; j < branchParticle->GetEntriesFast(); j++)
        {     
         
          gen = (GenParticle*) branchParticle->At(j);     
          candidate = factory->NewCandidate();

          candidate->Momentum = gen->P4();
          candidate->Position.SetXYZT(gen->X, gen->Y, gen->Z, gen->T*1.0E3*c_light);
          
          candidate->PID = gen->PID;
          candidate->Status = gen->Status;
     
          candidate->M1 = gen->M1;
          candidate->M2 = gen->M2;

          candidate->D1 = gen->D1;
          candidate->D2 = gen->D2;

          candidate->Charge = gen->Charge;
          candidate->Mass   = gen->Mass;
    
          allParticleOutputArray->Add(candidate);
 
          pdgCode = TMath::Abs(gen->PID);

          if(gen->Status == 1)
          {
            stableParticleOutputArray->Add(candidate);
          }
          else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15)
          {
            partonOutputArray->Add(candidate);
          }
        }
        
        modularDelphes->ProcessTask();

        treeWriter->Fill();

        modularDelphes->Clear();
        treeWriter->Clear();

        progressBar.Update(eventCounter, eventCounter);
        ++eventCounter;
      }
 
      progressBar.Update(eventCounter, eventCounter, kTRUE);
      progressBar.Finish();

      inputFile->Close();
    
      delete treeReader;
      
    }

    modularDelphes->FinishTask();
    treeWriter->Write();

    cout << "** Exiting..." << endl;

    delete modularDelphes;
    delete confReader;
    delete treeWriter;
    delete outputFile;
    delete chain;
    
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
