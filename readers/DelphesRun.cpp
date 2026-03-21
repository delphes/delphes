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
#include "classes/DelphesReader.h"
#include "modules/Delphes.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "DelphesRun";
  TStopwatch procStopWatch;
  Int_t i;

  if(argc < 3)
  {
    std::cout << " Usage: " << appName << " reader_type config_file output_file"
              << " [input_file(s)]" << std::endl;
    std::cout << " config_file - configuration file in Tcl format," << std::endl;
    std::cout << " output_file - output file in ROOT format," << std::endl;
    std::cout << " input_file(s) - input file(s) in HepMC format," << std::endl;
    std::cout << " with no input_file, or when input_file is -, read standard input." << std::endl;
    return 1;
  }

  signal(SIGINT, SignalHandler);

  gROOT->SetBatch();

  try
  {
    const auto confReader = std::make_unique<ExRootConfReader>();
    confReader->ReadFile(argv[2]);

    const auto reader = DelphesReaderFactory::Get().Build(argv[1]);
    reader->SetMaxEvents(confReader->GetInt("::MaxEvents", 0));
    reader->SetSkipEvents(confReader->GetInt("::SkipEvents", 0));

    const auto modularDelphes = std::make_unique<Delphes>("Delphes");
    modularDelphes->SetConfReader(confReader.get());
    modularDelphes->SetOutputFile(argv[3]);
    modularDelphes->SetReader(reader.get());

    modularDelphes->InitTask();

    i = 4;
    do
    {
      if(interrupted) break;

      if(i == argc || strncmp(argv[i], "-", 2) == 0)
      {
        std::cout << "** Reading standard input" << std::endl;
        throw;
        //auto *inputFile = stdin;
        //Long64_t length = -1;
      }
      else
      {
        std::cout << "** Reading " << argv[i] << std::endl;
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

    std::cout << "** Exiting..." << std::endl;

    return 0;
  }
  catch(const std::runtime_error &e)
  {
    std::cerr << "** ERROR: " << e.what() << std::endl;
    return 1;
  }
}
