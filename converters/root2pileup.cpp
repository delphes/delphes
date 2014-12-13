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
#include <string>

#include <signal.h>

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TClonesArray.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesPileUpWriter.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

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
  char appName[] = "root2pileup";
  stringstream message;
  TChain *inputChain = 0;
  ExRootTreeReader *treeReader = 0;
  TClonesArray *branchParticle = 0;
  TIterator *itParticle = 0;
  GenParticle *particle = 0;
  DelphesPileUpWriter *writer = 0;
  Long64_t entry, allEntries;
  Int_t i;

  if(argc < 3)
  {
    cout << " Usage: " << appName << " output_file" << " input_file(s)" << endl;
    cout << " output_file - output binary pile-up file," << endl;
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
    inputChain = new TChain("Delphes");
    for(i = 2; i < argc && !interrupted; ++i)
    {
      inputChain->Add(argv[i]);
    }

    treeReader = new ExRootTreeReader(inputChain);
    branchParticle = treeReader->UseBranch("Particle");
    itParticle = branchParticle->MakeIterator();

    writer = new DelphesPileUpWriter(argv[1]);

    allEntries = treeReader->GetEntries();
    cout << "** Input file(s) contain(s) " << allEntries << " events" << endl;

    if(allEntries > 0)
    {
      ExRootProgressBar progressBar(allEntries - 1);
      // Loop over all events in the input file
      for(entry = 0; entry < allEntries && !interrupted; ++entry)
      {
        if(!treeReader->ReadEntry(entry))
        {
          cerr << "** ERROR: cannot read event " << entry << endl;
          break;
        }

        itParticle->Reset();
        while((particle = static_cast<GenParticle*>(itParticle->Next())))
        {
          writer->WriteParticle(particle->PID,
            particle->X, particle->Y, particle->Z, particle->T,
            particle->Px, particle->Py, particle->Pz, particle->E);
        }
        
        writer->WriteEntry();

        progressBar.Update(entry);
      }
      progressBar.Finish();

      writer->WriteIndex();
    }

    cout << "** Exiting..." << endl;

    delete writer;
    delete itParticle;
    delete treeReader;
    delete inputChain;
    return 0;
  }
  catch(runtime_error &e)
  {
    if(writer) delete writer;
    if(itParticle) delete itParticle;
    if(treeReader) delete treeReader;
    if(inputChain) delete inputChain;
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}
