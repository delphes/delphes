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

#include <TApplication.h>
#include <TFile.h>
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

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  try
  {
    const auto outputFile = std::make_unique<TFile>(argv[3], "CREATE");

    if(!outputFile)
    {
      std::ostringstream message;
      message << "can't create output file " << argv[3];
      throw std::runtime_error(message.str());
    }

    const auto treeWriter = std::make_unique<ExRootTreeWriter>(outputFile.get(), "Delphes");

    ExRootTreeBranch *branchEvent = treeWriter->NewBranch("Event", HepMCEvent::Class()),
                     *branchWeight = treeWriter->NewBranch("Weight", Weight::Class());

    const auto confReader = std::make_unique<ExRootConfReader>();
    confReader->ReadFile(argv[2]);

    const auto modularDelphes = std::make_unique<Delphes>("Delphes");
    modularDelphes->SetConfReader(confReader.get());
    modularDelphes->SetTreeWriter(treeWriter.get());

    DelphesFactory *factory = modularDelphes->GetFactory();
    CandidatesCollection allParticleOutputArray = modularDelphes->ExportArray("allParticles"),
                         stableParticleOutputArray = modularDelphes->ExportArray("stableParticles"),
                         partonOutputArray = modularDelphes->ExportArray("partons");

    const auto reader = DelphesReaderFactory::Get().Build(argv[1]);
    reader->SetMaxEvents(confReader->GetInt("::MaxEvents", 0));
    reader->SetSkipEvents(confReader->GetInt("::SkipEvents", 0));

    modularDelphes->InitTask();

    i = 3;
    do
    {
      if(interrupted) break;

      if(i == argc || strncmp(argv[i], "-", 3) == 0)
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
      treeWriter->Clear();
      modularDelphes->Clear();
      reader->Clear();
      while(reader->ReadEvent(factory, allParticleOutputArray, stableParticleOutputArray, partonOutputArray) && !interrupted)
      {
        procStopWatch.Start();
        modularDelphes->ProcessTask();
        procStopWatch.Stop();

        reader->AnalyzeEvent(branchEvent, &procStopWatch);
        reader->AnalyzeWeight(branchWeight);

        treeWriter->Fill();

        treeWriter->Clear();

        modularDelphes->Clear();
        reader->Clear();
      }

      ++i;
    } while(i < argc);

    modularDelphes->FinishTask();
    treeWriter->Write();

    std::cout << "** Exiting..." << std::endl;

    return 0;
  }
  catch(const std::runtime_error &e)
  {
    std::cerr << "** ERROR: " << e.what() << std::endl;
    return 1;
  }
}
