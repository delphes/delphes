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


/*
########################################################################


This simple example shows how to use Delphes with an external FastJet installation.
Events in HepMC format are read via the DelphesHepMC reader. 

In order to run this example you first need to set the paths to your Delphes, FastJet
and ROOT installations (DELPHES_DIR, FASTJET_DIR and ROOT_DIR):

DELPHES_DIR=<path to Delphes installation>
FASTJET_DIR=<path to FastJet installation>
ROOT_DIR=<path to ROOT installation>

Then run the following commands to build the executable:

DELPHES_LIB="-Wl,-rpath,$DELPHES_DIR -L$DELPHES_DIR -lDelphesNoFastJet"

FASTJET_INC=`$FASTJET_DIR/bin/fastjet-config --cxxflags`
FASTJET_LIB=`$FASTJET_DIR/bin/fastjet-config --libs`

ROOT_INC=`$ROOT_DIR/bin/root-config --incdir`
ROOT_LIB=`$ROOT_DIR/bin/root-config --libs`

CXXFLAGS="$FASTJET_INC -I$ROOT_INC -I$DELPHES_DIR -I$DELPHES_DIR/external"
LDFLAGS="$FASTJET_LIB $ROOT_LIB -lEG $DELPHES_LIB"

g++ $CXXFLAGS examples/ExternalFastJet/ExternalFastJetHepMC.cpp $LDFLAGS -o ExternalFastJetHepMC

Then run (you need an event file in HepMC format):

./ExternalFastJetHepMC cards/delphes_card_CMS_NoFastJet.tcl file.hepmc


########################################################################
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

// #include "fastjet/contrib/Nsubjettiness.hh"
// #include "fastjet/contrib/Njettiness.hh"
// #include "fastjet/contrib/NjettinessPlugin.hh"
// #include "fastjet/contrib/WinnerTakeAllRecombiner.hh"

using namespace std;
using namespace fastjet;
// using namespace fastjet::contrib;

//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "ExternalFastJetHepMC";
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
//  JetDefinition::Recombiner *recomb = 0;
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

//    recomb = new WinnerTakeAllRecombiner();
//    definition = new JetDefinition(antikt_algorithm, 0.5, recomb, Best);
  
    definition = new JetDefinition(antikt_algorithm, 0.5);
    
    
    // Define your input candidates to fastjet (by default particle-flow objects).
    // If you want pure calorimeter towers change "EFlowMerger/eflow" into "Calorimeter/towers":
     
    inputArray = modularDelphes->ImportArray("EFlowMerger/eflow");

    inputIterator = inputArray->MakeIterator();


    // start reading hepmc file

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
	
	 // loop over events
	if(reader->EventReady())
        {
          ++eventCounter;

          if(eventCounter > skipEvents)
          {
            
	    // run delphes reconstruction
	    modularDelphes->ProcessTask();
            
            inputList.clear();
            inputIterator->Reset();
	
	
	    // pass delphes candidates to fastjet clustering  
             while((candidate = static_cast<Candidate*>(inputIterator->Next())))
            {
              momentum = candidate->Momentum;
              jet = PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
              inputList.push_back(jet);
            }
           
	    // run fastjet clustering 
	    ClusterSequence sequence(inputList, *definition);
            outputList.clear();
            outputList = sorted_by_pt(sequence.inclusive_jets(0.0));

            
	    // Prints for the first event:
            //  - the description of the algorithm used
            //  - show the inclusive jets as
            //      {index, rapidity, phi, pt}
            //----------------------------------------------------------
           
	    if(eventCounter == skipEvents + 1)
            {
              cout << "Ran " << definition->description() << endl;

              // label the columns
              printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");

              // print out the details for each jet
              for (unsigned int i = 0; i < outputList.size(); i++) {
                printf("%5u %15.8f %15.8f %15.8f\n",
                       i, outputList[i].rap(), outputList[i].phi(),
                       outputList[i].perp());
              }
            }
          }

          modularDelphes->Clear();
          reader->Clear();
        }
      } // end of event loop

      if(inputFile != stdin) fclose(inputFile);

      ++i;
    }
    while(i < argc);

    modularDelphes->FinishTask();

    cout << "** Exiting..." << endl;

    delete definition;
//    delete recomb;

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
