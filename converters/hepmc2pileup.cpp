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
#include <TLorentzVector.h>
#include <TROOT.h>

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesHepMC2Reader.h"
#include "classes/DelphesPileUpWriter.h"

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
  char appName[] = "hepmc2pileup";
  stringstream message;
  Int_t i;

  if(argc < 2)
  {
    cout << " Usage: " << appName << " output_file"
         << " [input_file(s)]" << endl;
    cout << " output_file - output binary pile-up file," << endl;
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
    const auto writer = std::make_unique<DelphesPileUpWriter>(argv[1]);

    const auto factory = std::make_unique<DelphesFactory>();

    const auto reader = std::make_unique<DelphesHepMC2Reader>();
    reader->SetFactory(factory.get());

    i = 2;
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
      factory->Clear();
      reader->Clear();
      CandidatesCollection stableParticleOutputArray = factory->Attach<std::vector<Candidate *> >("stableParticles");
      while(reader->ReadEvent() && !interrupted)
      {
        for(const auto &candidate : *stableParticleOutputArray)
        {
          const TLorentzVector &position = candidate->Position;
          const TLorentzVector &momentum = candidate->Momentum;
          writer->WriteParticle(candidate->PID,
            position.X(), position.Y(), position.Z(), position.T(),
            momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
        }

        writer->WriteEntry();

        factory->Clear();
        reader->Clear();
      }
      ++i;
    } while(i < argc);

    writer->WriteIndex();

    cout << "** Exiting..." << endl;

    return 0;
  }
  catch(runtime_error &e)
  {
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}
