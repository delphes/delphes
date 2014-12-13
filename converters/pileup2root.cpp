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

#include <stdlib.h>
#include <signal.h>
#include <stdio.h>

#include <rpc/types.h>
#include <rpc/xdr.h>

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TObjArray.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

#include "classes/DelphesStream.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesPileUpReader.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

using namespace std;

//------------------------------------------------------------------------------

void ProcessEvent(DelphesPileUpReader *reader, ExRootTreeBranch *branch)
{
  GenParticle *particle;
  Int_t pid;
  Float_t x, y, z, t;
  Float_t px, py, pz, e;
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  TParticlePDG *pdgParticle;
  TLorentzVector momentum;
  Double_t pt, signPz, cosTheta, eta, rapidity;

  while(reader->ReadParticle(pid, x, y, z, t, px, py, pz, e))
  {
    particle = static_cast<GenParticle*>(branch->NewEntry());

    particle->PID = pid;
    particle->X = x;
    particle->Y = y;
    particle->Z = z;
    particle->T = t;
    particle->Px = px;
    particle->Py = py;
    particle->Pz = pz;
    particle->E = e;

    particle->Status = 1;
    particle->IsPU = 1;

    particle->M1 = -1;
    particle->M2 = -1;

    particle->D1 = -1;
    particle->D2 = -1;

    pdgParticle = pdg->GetParticle(pid);
    particle->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;

    particle->Mass = pdgParticle ? pdgParticle->Mass() : -999.9;

    momentum.SetPxPyPzE(px, py, pz, e);
    pt = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signPz*999.9 : momentum.Eta());
    rapidity = (cosTheta == 1.0 ? signPz*999.9 : momentum.Rapidity());

    particle->Eta = eta;
    particle->Phi = momentum.Phi();
    particle->PT = pt;

    particle->Rapidity = rapidity;
  }
}

//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "pileup2root";
  stringstream message;
  TFile *outputFile = 0;
  ExRootTreeWriter *treeWriter = 0;
  ExRootTreeBranch *branchParticle = 0;
  DelphesPileUpReader *reader = 0;
  Long64_t entry, allEntries;

  if(argc != 3)
  {
    cout << " Usage: " << appName << " output_file" << " input_file" << endl;
    cout << " output_file - output file in ROOT format," << endl;
    cout << " input_file - input binary pile-up file." << endl;
    return 1;
  }

  signal(SIGINT, SignalHandler);

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  try
  {
    outputFile = TFile::Open(argv[1], "CREATE");

    if(outputFile == NULL)
    {
      message << "can't open " << argv[1];
      throw runtime_error(message.str());
    }

    cout << "** Reading " << argv[2] << endl;

    reader = new DelphesPileUpReader(argv[2]);
    allEntries = reader->GetEntries();

    cout << "** Input file contains " << allEntries << " events" << endl;

    if(allEntries > 0)
    {
      treeWriter = new ExRootTreeWriter(outputFile, "Delphes");
      branchParticle = treeWriter->NewBranch("Particle", GenParticle::Class());

      ExRootProgressBar progressBar(allEntries - 1);
      // Loop over all events
      for(entry = 0; entry < allEntries && !interrupted; ++entry)
      {
        if(!reader->ReadEntry(entry))
        {
          cerr << "** ERROR: cannot read event " << entry << endl;
          break;
        }

        treeWriter->Clear();
        ProcessEvent(reader, branchParticle);
        treeWriter->Fill();

        progressBar.Update(entry);
      }
      treeWriter->Write();

      progressBar.Finish();

      delete treeWriter;
    }
    delete reader;

    cout << "** Exiting..." << endl;

    return 0;
  }
  catch(runtime_error &e)
  {
    if(treeWriter) delete treeWriter;
    if(reader) delete reader;
    if(outputFile) delete outputFile;
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}


