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

#include <signal.h>

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesReader.h"
#include "classes/DelphesTCLConfReader.h"
#include "modules/Delphes.h"

#include <ExRootAnalysis/ExRootProgressBar.h>

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

  if(argc != 4)
  {
    std::cout << " Usage: " << appName << " config_file"
              << " pythia_card"
              << " output_file" << std::endl
              << " config_file - configuration file in Tcl format," << std::endl
              << " pythia_card - Pythia8 configuration file," << std::endl
              << " output_file - output file in ROOT format." << std::endl;
    return 1;
  }

  signal(SIGINT, SignalHandler);

  try
  {
    const auto confReader = std::make_unique<DelphesTCLConfReader>();
    confReader->ReadFile(argv[1]);

    const auto modularDelphes = std::make_unique<Delphes>("Delphes");
    modularDelphes->SetConfReader(confReader.get());
    modularDelphes->SetOutputFile(argv[3]);

    const auto reader = DelphesReaderFactory::Get().Build("Pythia8", DelphesParameters{}.Set("pythiaConfigFile", argv[2]));
    const auto userParams = confReader->Parameters();
    reader->SetMaxEvents(userParams.Get<int>("MaxEvents", 0));
    reader->SetSkipEvents(userParams.Get<int>("SkipEvents", 0));

    modularDelphes->InitTask();

    // Loop over all objects
    modularDelphes->Clear();
    reader->Clear();
    while(reader->ReadEvent() && !interrupted)
    {
      modularDelphes->ProcessTask();
      modularDelphes->Clear();
      reader->Clear();
    }
    // Loop over all events

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
