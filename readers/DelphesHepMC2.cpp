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

#include <signal.h>

#include <TROOT.h>
#include <TStopwatch.h>

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesHepMC2Reader.h"
#include "modules/Delphes.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

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
  char appName[] = "DelphesHepMC2";
  stringstream message;
  TStopwatch procStopWatch;
  Int_t i;

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

  try
  {
    const auto confReader = std::make_unique<ExRootConfReader>();
    confReader->ReadFile(argv[1]);

    const auto modularDelphes = std::make_unique<Delphes>("Delphes");
    modularDelphes->SetConfReader(confReader.get());

    const auto reader = std::make_unique<DelphesHepMC2Reader>();
    reader->SetMaxEvents(confReader->GetInt("::MaxEvents", 0));
    reader->SetSkipEvents(confReader->GetInt("::SkipEvents", 0));

    modularDelphes->InitTask();

    i = 3;
    do
    {
      if(interrupted) break;

      if(i == argc || strncmp(argv[i], "-", 2) == 0)
      {
        cout << "** Reading standard input" << endl;
        //inputFile = stdin;
        //length = -1;
        throw;
      }
      else
      {
        cout << "** Reading " << argv[i] << endl;
        reader->LoadInputFile(argv[i]);
      }

      // Loop over all objects
      modularDelphes->Clear();
      reader->Clear();
      while(reader->ReadEvent() && !interrupted)
      {
        procStopWatch.Start();
        modularDelphes->ProcessTask();
        procStopWatch.Stop();

        reader->AnalyzeEvent(&procStopWatch);

        modularDelphes->Clear();
        reader->Clear();
      }
      ++i;
    } while(i < argc);

    modularDelphes->FinishTask();

    cout << "** Exiting..." << endl;

    return 0;
  }
  catch(runtime_error &e)
  {
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}
