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
#include <sstream>
#include <stdexcept>

#include <csignal>

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
#include "classes/DelphesLHEFReader.h"
#include "classes/DelphesPythia8Reader.h"
#include "modules/Delphes.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "classes/DelphesTreeWriter.h"

using namespace std;

//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "DelphesPythia8";
  stringstream message;
  FILE *inputFileLHEF = 0;
  TFile *outputFile = 0;
  TStopwatch readStopWatch, procStopWatch;
  DelphesTreeWriter *treeWriter = 0;
  ExRootTreeBranch *branchEvent = 0, *branchWeight = 0;
  ExRootTreeBranch *branchEventLHEF = 0, *branchWeightLHEF = 0;
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  DelphesFactory *factory = 0;
  TObjArray *stableParticleOutputArray = 0, *allParticleOutputArray = 0, *partonOutputArray = 0;
  TObjArray *stableParticleOutputArrayLHEF = 0, *allParticleOutputArrayLHEF = 0, *partonOutputArrayLHEF = 0;
  DelphesPythia8Reader *reader = 0;
  DelphesLHEFReader *readerLHEF = 0;
  Long64_t eventCounter;
  Bool_t readLHEF;

  if(argc != 4)
  {
    cout << " Usage: " << appName << " config_file pythia_card output_file" << endl;
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

    treeWriter = new DelphesTreeWriter(outputFile, "Delphes");

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

    reader = new DelphesPythia8Reader;
    reader->OpenInputFile(argv[2]);

    readLHEF = reader->ReadLHEF();
    if(readLHEF)
    {
      inputFileLHEF = fopen(reader->FileNameLHEF().c_str(), "r");
      if(inputFileLHEF == NULL)
      {
        message << "can't open LHEF file " << reader->FileNameLHEF();
        throw runtime_error(message.str());
      }

      readerLHEF = new DelphesLHEFReader;
      readerLHEF->SetInputFile(inputFileLHEF);

      branchEventLHEF = treeWriter->NewBranch("EventLHEF", LHEFEvent::Class());
      branchWeightLHEF = treeWriter->NewBranch("WeightLHEF", LHEFWeight::Class());

      allParticleOutputArrayLHEF = modularDelphes->ExportArray("allParticlesLHEF");
      stableParticleOutputArrayLHEF = modularDelphes->ExportArray("stableParticlesLHEF");
      partonOutputArrayLHEF = modularDelphes->ExportArray("partonsLHEF");
    }

    modularDelphes->Init();


    ExRootProgressBar progressBar(-1);

    // Loop over all events
    eventCounter = 0;
    treeWriter->Clear();
    modularDelphes->Clear();
    reader->Clear();
    if(readLHEF) readerLHEF->Clear();
    readStopWatch.Start();
    while(reader->ReadEvent(factory, allParticleOutputArray, stableParticleOutputArray, partonOutputArray) && !interrupted)
    {
      ++eventCounter;

      if(readLHEF)
      {
        while(eventCounter < reader->EventNumber())
        {
          ++eventCounter;
          readerLHEF->SkipEvent();
        }
        readerLHEF->ReadEvent(factory, allParticleOutputArrayLHEF, stableParticleOutputArrayLHEF, partonOutputArrayLHEF);
      }

      readStopWatch.Stop();

      if(reader->EventReady())
      {
        procStopWatch.Start();

        modularDelphes->Process();
        procStopWatch.Stop();

        reader->AnalyzeEvent(branchEvent, eventCounter, &readStopWatch, &procStopWatch);
        reader->AnalyzeWeight(branchWeight);

        if(readLHEF)
        {
          readerLHEF->AnalyzeEvent(branchEventLHEF, eventCounter, &readStopWatch, &procStopWatch);
          readerLHEF->AnalyzeWeight(branchWeightLHEF);
        }

        treeWriter->Fill();
      }

      treeWriter->Clear();
      modularDelphes->Clear();
      reader->Clear();
      if(readLHEF) readerLHEF->Clear();

      readStopWatch.Start();
      progressBar.Update(eventCounter, eventCounter);
    }

    progressBar.Update(eventCounter, eventCounter, kTRUE);
    progressBar.Finish();

    modularDelphes->Finish();
    treeWriter->Write();

    cout << "** Exiting..." << endl;

    delete readerLHEF;
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
