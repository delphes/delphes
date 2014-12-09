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

#include <vector>

#include <signal.h>

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

#include "modules/Delphes.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesHepMCReader.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "fastjet/plugins/SISCone/fastjet/SISConePlugin.hh"
#include "fastjet/plugins/CDFCones/fastjet/CDFMidPointPlugin.hh"
#include "fastjet/plugins/CDFCones/fastjet/CDFJetCluPlugin.hh"

#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"
#include "fastjet/contribs/Nsubjettiness/Njettiness.hh"
#include "fastjet/contribs/Nsubjettiness/NjettinessPlugin.hh"
#include "fastjet/contribs/Nsubjettiness/WinnerTakeAllRecombiner.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "StandaloneHepMC";
  stringstream message;
  FILE *inputFile = 0;
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  DelphesFactory *factory = 0;
  TObjArray *stableParticleOutputArray = 0, *allParticleOutputArray = 0, *partonOutputArray = 0;
  DelphesHepMCReader *reader = 0;
  Int_t i, maxEvents, skipEvents;
  Long64_t length, eventCounter;
  
  TObjArray *inputArray = 0;
  TIterator *inputIterator = 0;
  Candidate *candidate = 0;
  TLorentzVector momentum;

  JetDefinition *definition = 0;
  JetDefinition::Recombiner *recomb = 0;
  vector<PseudoJet> inputList, outputList;
  PseudoJet jet;

  if(argc < 2)
  {
    cout << " Usage: " << appName << " config_file" << " [input_file(s)]" << endl;
    cout << " config_file - configuration file in Tcl format," << endl;
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

    factory = modularDelphes->GetFactory();
    allParticleOutputArray = modularDelphes->ExportArray("allParticles");
    stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
    partonOutputArray = modularDelphes->ExportArray("partons");

    reader = new DelphesHepMCReader;

    modularDelphes->InitTask();

    ClusterSequence::print_banner();
    recomb = new WinnerTakeAllRecombiner();
    definition = new JetDefinition(antikt_algorithm, 0.5, recomb, Best);

    inputArray = modularDelphes->ImportArray("Calorimeter/towers");
    inputIterator = inputArray->MakeIterator();

    i = 2;
    do
    {
      if(interrupted) break;

      if(i == argc || strncmp(argv[i], "-", 2) == 0)
      {
        cout << "** Reading standard input" << endl;
        inputFile = stdin;
        length = -1;
      }
      else
      {
        cout << "** Reading " << argv[i] << endl;
        inputFile = fopen(argv[i], "r");

        if(inputFile == NULL)
        {
          message << "can't open " << argv[i];
          throw runtime_error(message.str());
        }

        fseek(inputFile, 0L, SEEK_END);
        length = ftello(inputFile);
        fseek(inputFile, 0L, SEEK_SET);

        if(length <= 0)
        {
          fclose(inputFile);
          ++i;
          continue;
        }
      }

      reader->SetInputFile(inputFile);

      // Loop over all objects
      eventCounter = 0;
      modularDelphes->Clear();
      reader->Clear();
      while((maxEvents <= 0 || eventCounter - skipEvents < maxEvents) &&
        reader->ReadBlock(factory, allParticleOutputArray,
        stableParticleOutputArray, partonOutputArray) && !interrupted)
      {
        if(reader->EventReady())
        {
          ++eventCounter;

          if(eventCounter > skipEvents)
          {
            modularDelphes->ProcessTask();
            
            inputList.clear();
            inputIterator->Reset();
            while((candidate = static_cast<Candidate*>(inputIterator->Next())))
            {
              momentum = candidate->Momentum;
              jet = PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
              inputList.push_back(jet);
            }
            ClusterSequence sequence(inputList, *definition);
            outputList.clear();
            outputList = sorted_by_pt(sequence.inclusive_jets(10.0));
          }

          modularDelphes->Clear();
          reader->Clear();
        }
      }

      if(inputFile != stdin) fclose(inputFile);

      ++i;
    }
    while(i < argc);

    modularDelphes->FinishTask();

    cout << "** Exiting..." << endl;

    delete reader;
    delete modularDelphes;
    delete confReader;

    return 0;
  }
  catch(runtime_error &e)
  {
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}
