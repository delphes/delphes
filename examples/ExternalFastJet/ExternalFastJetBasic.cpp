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
Events are specified via the multidimensional array "EVENTS" (for an example reading 
an HepMC file see ExternalFastJetHepMC.cpp). 

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

g++ $CXXFLAGS examples/ExternalFastJet/ExternalFastJetBasic.cpp $LDFLAGS -o ExternalFastJetBasic

Then run:

./ExternalFastJetBasic cards/delphes_card_CMS_NoFastJet.tcl


########################################################################
*/

#include <stdexcept>
#include <iostream>
#include <sstream>

#include <stdio.h>

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

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------------------

const int NEVENTS = 2;
const int NPARTICLES = 4;

float EVENTS[NEVENTS][NPARTICLES][11] =
{
  {
    {-211, 1, -2.91e-01,  1.99e-01, -1.30e+00, 1.35e+00, 1.3957e-01, 0.0, 0.0, 0.0, 0.0},
    { 211, 1, -2.62e-02,  1.90e-01,  1.42e+00, 1.44e+00, 1.3957e-01, 0.0, 0.0, 0.0, 0.0},
    {-211, 1, -3.08e-02, -2.30e-01,  5.13e+00, 5.14e+00, 1.3957e-01, 0.0, 0.0, 0.0, 0.0},
    { 211, 1,  1.21e-01, -3.52e-01,  2.12e+01, 2.12e+01, 1.3957e-01, 0.0, 0.0, 0.0, 0.0}
  },
  {
    {22, 1, -7.80e-02, -9.70e-02, 2.03e+01, 2.03e+01, 0.0, 0.0, 0.0, 0.0, 0.0},
    {22, 1, -1.16e-02, -2.08e-02, 4.89e-01, 4.90e-01, 0.0, 0.0, 0.0, 0.0, 0.0},
    {22, 1,  5.64e-02, -1.34e-02, 5.81e+00, 5.82e+00, 0.0, 0.0, 0.0, 0.0, 0.0},
    {22, 1,  2.35e-01, -1.04e-01, 9.12e+00, 9.12e+00, 0.0, 0.0, 0.0, 0.0, 0.0}
  }
};

//---------------------------------------------------------------------------


// this function converts input event array into Delphes candidates (defined below)

void ConvertInput(Int_t event, DelphesFactory *factory, TObjArray *allParticleOutputArray, TObjArray *stableParticleOutputArray, TObjArray *partonOutputArray);


//----------------------------------------------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  // Declaration of variables
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  DelphesFactory *factory = 0;
  TObjArray *allParticleOutputArray = 0;
  TObjArray *stableParticleOutputArray = 0;
  TObjArray *partonOutputArray = 0;

  Int_t event;

  TObjArray *inputArray = 0;
  TIterator *inputIterator = 0;
  Candidate *candidate = 0;
  TLorentzVector momentum;

  JetDefinition *definition = 0;
  vector<PseudoJet> inputList, outputList;
  PseudoJet jet;

  gROOT->SetBatch();

  int appargc = 1;
  char appName[] = "ExternalFastJetBasic";
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  if(argc != 2)
  {
    cout << " Usage: " << appName << " config_file" << endl;
    cout << " config_file - configuration file in Tcl format." << endl;
    return 1;
  }

  try
  {
    // Initialization
    confReader = new ExRootConfReader;
    confReader->ReadFile(argv[1]);

    modularDelphes = new Delphes("Delphes");
    modularDelphes->SetConfReader(confReader);

    factory = modularDelphes->GetFactory();

    allParticleOutputArray = modularDelphes->ExportArray("allParticles");
    stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
    partonOutputArray = modularDelphes->ExportArray("partons");

    modularDelphes->InitTask();

    
    // fastjet definition
    ClusterSequence::print_banner();
    definition = new JetDefinition(antikt_algorithm, 0.5);
    
    // Define your input candidates to fastjet (by default particle-flow objects).
    // If you want pure calorimeter towers change "EFlowMerger/eflow" into "Calorimeter/towers":
     
    inputArray = modularDelphes->ImportArray("EFlowMerger/eflow");
      
    inputIterator = inputArray->MakeIterator();

    // Event loop
    for(event = 0; event < NEVENTS; ++event)
    {
      modularDelphes->Clear();
      
      // convert EVENT input array into Delphes internal format
      ConvertInput(event, factory, allParticleOutputArray, stableParticleOutputArray, partonOutputArray);
      
      // run Delphes reconstruction
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
     
      // run clustering 
      ClusterSequence sequence(inputList, *definition);
      outputList.clear();
      outputList = sorted_by_pt(sequence.inclusive_jets(0.0));

      // tell the user what was done
      //  - the description of the algorithm used
      //  - show the inclusive jets as
      //      {index, rapidity, phi, pt}
      //----------------------------------------------------------
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

    // Finalization
    modularDelphes->FinishTask();
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


// ------------------------------------------------------------------------------------------------------------------------------------

void ConvertInput(Int_t event, DelphesFactory *factory, TObjArray *allParticleOutputArray, TObjArray *stableParticleOutputArray, TObjArray *partonOutputArray)
{
  Int_t particle;

  Candidate *candidate;
  TDatabasePDG *pdg;
  TParticlePDG *pdgParticle;
  Int_t pdgCode;

  Int_t pid, status;
  Double_t px, py, pz, e, mass;
  Double_t x, y, z, t;

  pdg = TDatabasePDG::Instance();

  for(particle = 0; particle < NPARTICLES; ++particle)
  {
    pid = EVENTS[event][particle][0];
    status = EVENTS[event][particle][1];
    px = EVENTS[event][particle][2];
    py = EVENTS[event][particle][3];
    pz = EVENTS[event][particle][4];
    e = EVENTS[event][particle][5];
    mass = EVENTS[event][particle][6];
    x = EVENTS[event][particle][7];
    y = EVENTS[event][particle][8];
    z = EVENTS[event][particle][9];
    t = EVENTS[event][particle][10];

    candidate = factory->NewCandidate();

    candidate->PID = pid;
    pdgCode = TMath::Abs(candidate->PID);

    candidate->Status = status;

    pdgParticle = pdg->GetParticle(pid);
    candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
    candidate->Mass = mass;

    candidate->Momentum.SetPxPyPzE(px, py, pz, e);

    candidate->Position.SetXYZT(x, y, z, t);

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

