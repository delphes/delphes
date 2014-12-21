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
#include <fstream>
#include <sstream>
#include <string>

#include <stdlib.h>
#include <signal.h>
#include <stdio.h>

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TClonesArray.h"

#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

using namespace std;

/*
LHC Olympics format discription from http://www.jthaler.net/olympicswiki/doku.php?id=lhc_olympics:data_file_format

    * The first column of each row is just a counter that labels the object.
    * The event begins with a row labelled "0"; this row contains the event number and the triggering information. The last row of the event is always the missing transverse momentum (MET).
    * The second column of each row gives the type of object being listed [0, 1, 2, 3, 4, 6 = photon, electron, muon, hadronically-decaying tau, jet, missing transverse energy].
    * The next three columns give the pseudorapidity, the azimuthal angle, and the transverse momentum of the object.
    * The sixth column gives the invariant mass of the object.
    * The seventh column gives the number of tracks associated with the object; in the case of a lepton, this number is multiplied by the charge of the lepton.
    * The eighth column is 1 or 2 for a jet that has been "tagged" as containing a b-quark (actually a heavy flavor tag that sometimes indicates c-quarks), otherwise it is 0. For muons, the integer part of this number is the identity of the jet (see column 1) that is closest ot this muon in Delta R.
    * The ninth column is the ratio of the hadronic versus electromagnetic energy deposited in the calorimeter cells associated with the object. For muons to the left of the decimal point is the summed pT in a R=0.4 cone (excluding the muon). To the right of the decimal point is etrat, which is a percentage between .00 and .99. It is the ratio of the transverse energy in a 3x3 grid surrounding the muon to the pT of the muon.
*/

class LHCOWriter
{
public:
  LHCOWriter(ExRootTreeReader *treeReader, FILE *outputFile);
  ~LHCOWriter();

  void ProcessEvent();

private:

  void Reset();
  void Write();

  void AnalyseEvent();

  void AnalysePhotons();
  void AnalyseElectrons();
  void AnalyseMuons();
  void AnalyseTauJets();
  void AnalyseJets();

  void AnalyseMissingET();

  enum {kIntParamSize = 2, kDblParamSize = 9};
  Int_t fIntParam[kIntParamSize];
  Double_t fDblParam[kDblParamSize];

  Long64_t fTriggerWord, fEventNumber;

  ExRootTreeReader *fTreeReader;
  FILE *fOutputFile;

  TClonesArray *fBranchEvent;

  TClonesArray *fBranchTrack;
  TClonesArray *fBranchTower;

  TClonesArray *fBranchPhoton;
  TClonesArray *fBranchElectron;
  TClonesArray *fBranchMuon;
  TClonesArray *fBranchJet;
  TClonesArray *fBranchMissingET;

  TIterator *fItTrack;
  TIterator *fItTower;

  TIterator *fItPhoton;
  TIterator *fItElectron;
  TIterator *fItMuon;
  TIterator *fItJet;
};

//------------------------------------------------------------------------------

LHCOWriter::LHCOWriter(ExRootTreeReader *treeReader, FILE *outputFile) :
  fTriggerWord(0), fEventNumber(1), fTreeReader(0), fOutputFile(0),
  fBranchEvent(0), fBranchTrack(0), fBranchTower(0), fBranchPhoton(0),
  fBranchElectron(0), fBranchMuon(0), fBranchJet(0), fBranchMissingET(0)
{
  fTreeReader = treeReader;
  fOutputFile = outputFile;

  // information about reconstructed event
  fBranchEvent = fTreeReader->UseBranch("Event");
  // reconstructed tracks
  fBranchTrack = fTreeReader->UseBranch("Track");
  // calorimeter towers
  fBranchTower = fTreeReader->UseBranch("Tower");
  // reconstructed photons
  fBranchPhoton = fTreeReader->UseBranch("Photon");
  // reconstructed electrons
  fBranchElectron = fTreeReader->UseBranch("Electron");
  // reconstructed muons
  fBranchMuon = fTreeReader->UseBranch("Muon");
  // reconstructed jets
  fBranchJet = fTreeReader->UseBranch("Jet");
  // missing transverse energy
  fBranchMissingET = fTreeReader->UseBranch("MissingET");

  if(!fBranchEvent || !fBranchTrack || !fBranchTower || !fBranchPhoton ||
     !fBranchElectron || !fBranchMuon || !fBranchJet || !fBranchMissingET)
  {
    throw runtime_error("ROOT file doesn't contain all required branches");
  }

  fItTrack = fBranchTrack->MakeIterator();
  fItTower = fBranchTower->MakeIterator();
  fItPhoton = fBranchPhoton->MakeIterator();
  fItElectron = fBranchElectron->MakeIterator();
  fItMuon = fBranchMuon->MakeIterator();
  fItJet = fBranchJet->MakeIterator();
}

//------------------------------------------------------------------------------

LHCOWriter::~LHCOWriter()
{
}

//---------------------------------------------------------------------------

void LHCOWriter::ProcessEvent()
{
  fIntParam[0] = 0;

  AnalyseEvent();

  AnalysePhotons();
  AnalyseElectrons();
  AnalyseMuons();
  AnalyseTauJets();
  AnalyseJets();

  AnalyseMissingET();
}

//---------------------------------------------------------------------------

void LHCOWriter::Reset()
{
  int i;
  for(i = 1; i < kIntParamSize; ++i)
  {
    fIntParam[i] = 0;
  }

  for(i = 0; i < kDblParamSize; ++i)
  {
    fDblParam[i] = 0.0;
  }
}

//---------------------------------------------------------------------------

void LHCOWriter::Write()
{
  fprintf(fOutputFile, "%4d %4d %8.3f %8.3f %7.2f %7.2f %6.1f %6.1f %7.2f %6.1f %6.1f\n",
    fIntParam[0], fIntParam[1], fDblParam[0], fDblParam[1], fDblParam[2],
    fDblParam[3], fDblParam[4], fDblParam[5], fDblParam[6], fDblParam[7], fDblParam[8]);

  ++fIntParam[0];
}

//---------------------------------------------------------------------------

void LHCOWriter::AnalyseEvent()
{
  Event *element;

  element = static_cast<Event*>(fBranchEvent->At(0));

  fprintf(fOutputFile, "%4d %13lld %8d\n", 0, element->Number, 0);

  ++fIntParam[0];
}

//---------------------------------------------------------------------------

void LHCOWriter::AnalysePhotons()
{
  Photon *element;

  fItPhoton->Reset();
  while((element = static_cast<Photon*>(fItPhoton->Next())))
  {
    Reset();

    fIntParam[1] = 0;

    fDblParam[0] = element->Eta;
    fDblParam[1] = element->Phi;
    fDblParam[2] = element->PT;

    fDblParam[6] = element->EhadOverEem;

    Write();
  }
}

//---------------------------------------------------------------------------

void LHCOWriter::AnalyseElectrons()
{
  Electron *element;

  fItElectron->Reset();
  while((element = static_cast<Electron*>(fItElectron->Next())))
  {
    Reset();

    fIntParam[1] = 1;

    fDblParam[0] = element->Eta;
    fDblParam[1] = element->Phi;
    fDblParam[2] = element->PT;

    fDblParam[4] = element->Charge;

    fDblParam[6] = element->EhadOverEem;

    Write();
  }
}

//---------------------------------------------------------------------------

void LHCOWriter::AnalyseMuons()
{
  Muon *element;
  Track *track;
  Tower *tower;
  Jet *jet;
  Int_t muonCounter, tauCounter, jetCounter, minIndex;
  Float_t sumPT, sumET, ratET, jetDR, minDR;

  muonCounter = 0;
  fItMuon->Reset();
  while((element = static_cast<Muon*>(fItMuon->Next())))
  {
    Reset();

    sumPT = 0.0;
    fItTrack->Reset();
    while((track = static_cast<Track*>(fItTrack->Next())))
    {
      if(element->P4().DeltaR(track->P4()) < 0.5) sumPT += track->PT;
    }

    sumET = 0.0;
    fItTower->Reset();
    while((tower = static_cast<Tower*>(fItTower->Next())))
    {
      if(element->P4().DeltaR(tower->P4()) < 0.5) sumET += tower->ET;
    }

    tauCounter = 0;
    jetCounter = 0;
    minIndex = -1;
    minDR = 1.0E9;
    fItJet->Reset();
    while((jet = static_cast<Jet*>(fItJet->Next())))
    {
      if(jet->TauTag != 0)
      {
        ++tauCounter;
        continue;
      }

      jetDR = element->P4().DeltaR(jet->P4());
      if(jetDR < minDR)
      {
        minIndex = jetCounter;
        minDR = jetDR;
      }
      ++jetCounter;
    }

    fIntParam[1] = 2;

    fDblParam[0] = element->Eta;
    fDblParam[1] = element->Phi;
    fDblParam[2] = element->PT;

    fDblParam[3] = 0.11;

    fDblParam[4] = element->Charge;

    if(minIndex >= 0)
    {
      fDblParam[5] = fIntParam[0] + fBranchMuon->GetEntriesFast() - muonCounter + tauCounter + minIndex;
    }

    ratET = sumET/element->PT;
    fDblParam[6] = Float_t(TMath::Nint(sumPT)) + (ratET < 1.0 ? ratET : 0.99);

    Write();
    ++muonCounter;
  }
}

//---------------------------------------------------------------------------

void LHCOWriter::AnalyseTauJets()
{
  Jet *element;
  Track *track;
  Int_t counter;

  fItJet->Reset();
  while((element = static_cast<Jet*>(fItJet->Next())))
  {
    if(element->TauTag == 0) continue;

    Reset();

    counter = 0;
    fItTrack->Reset();
    while((track = static_cast<Track*>(fItTrack->Next())))
    {
      if(element->P4().DeltaR(track->P4()) < 0.5) ++counter;
    }

    fIntParam[1] = 3;

    fDblParam[0] = element->Eta;
    fDblParam[1] = element->Phi;
    fDblParam[2] = element->PT;
    fDblParam[3] = element->Mass;
    fDblParam[4] = counter * element->Charge;

    fDblParam[6] = element->EhadOverEem;

    Write();
  }
}

//---------------------------------------------------------------------------

void LHCOWriter::AnalyseJets()
{
  Jet *element;
  Track *track;
  Int_t counter;

  fItJet->Reset();
  while((element = static_cast<Jet*>(fItJet->Next())))
  {
    if(element->TauTag != 0) continue;

    Reset();

    counter = 0;
    fItTrack->Reset();
    while((track = static_cast<Track*>(fItTrack->Next())))
    {
      if(element->P4().DeltaR(track->P4()) < 0.5) ++counter;
    }

    fIntParam[1] = 4;

    fDblParam[0] = element->Eta;
    fDblParam[1] = element->Phi;
    fDblParam[2] = element->PT;
    fDblParam[3] = element->Mass;
    fDblParam[4] = counter;
    fDblParam[5] = element->BTag;
    fDblParam[6] = element->EhadOverEem;

    Write();
  }
}

//---------------------------------------------------------------------------

void LHCOWriter::AnalyseMissingET()
{
  MissingET *element;

  element = static_cast<MissingET*>(fBranchMissingET->At(0));

  Reset();

  fIntParam[1] = 6;

  fDblParam[1] = element->Phi;
  fDblParam[2] = element->MET;

  Write();
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
  char appName[] = "root2lhco";
  stringstream message;
  FILE *outputFile = 0;
  TChain *inputChain = 0;
  LHCOWriter *writer = 0;
  ExRootTreeReader *treeReader = 0;
  Long64_t entry, allEntries;

  if(argc < 2 || argc > 3)
  {
    cerr << " Usage: " << appName << " input_file" << " [output_file]" << endl;
    cerr << " input_file - input file in ROOT format," << endl;
    cerr << " output_file - output file in LHCO format," << endl;
    cerr << " with no output_file, or when output_file is -, write to standard output." << endl;
    return 1;
  }

  signal(SIGINT, SignalHandler);

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  try
  {
    cerr << "** Reading " << argv[1] << endl;
    inputChain = new TChain("Delphes");
    inputChain->Add(argv[1]);

    ExRootTreeReader *treeReader = new ExRootTreeReader(inputChain);

    if(argc == 2 || strcmp(argv[2], "-") == 0)
    {
      outputFile = stdout;
    }
    else
    {
      outputFile = fopen(argv[2], "w");

      if(outputFile == NULL)
      {
        message << "can't open " << argv[2];
        throw runtime_error(message.str());
      }
    }

    fprintf(outputFile, "   #  typ      eta      phi      pt    jmas   ntrk   btag  had/em   dum1   dum2\n");

    allEntries = treeReader->GetEntries();
    cerr << "** Input file contains " << allEntries << " events" << endl;

    if(allEntries > 0)
    {
      // Create LHC Olympics converter:
      writer = new LHCOWriter(treeReader, outputFile);

      ExRootProgressBar progressBar(allEntries - 1);
      // Loop over all events
      for(entry = 0; entry < allEntries && !interrupted; ++entry)
      {
        if(!treeReader->ReadEntry(entry))
        {
          cerr << "** ERROR: cannot read event " << entry << endl;
          break;
        }

        writer->ProcessEvent();

        progressBar.Update(entry);
      }
      progressBar.Finish();

      delete writer;
    }

    cerr << "** Exiting..." << endl;

    if(outputFile != stdout) fclose(outputFile);
    delete treeReader;
    delete inputChain;
    return 0;
  }
  catch(runtime_error &e)
  {
    if(writer) delete writer;
    if(treeReader) delete treeReader;
    if(inputChain) delete inputChain;
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}


