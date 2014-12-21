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

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TObjArray.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

#include "modules/Delphes.h"
#include "classes/DelphesStream.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

using namespace std;

static const int kBufferSize  = 1024;

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

//------------------------------------------------------------------------------

class LHCOConverter
{
public:
  LHCOConverter(TFile *outputFile);
  ~LHCOConverter();

  void Write();

  Bool_t ReadLine(FILE *inputFile);

private:

  void AddMissingEvents();

  void AnalyseEvent(ExRootTreeBranch *branch);

  void AnalysePhoton(ExRootTreeBranch *branch);
  void AnalyseElectron(ExRootTreeBranch *branch);
  void AnalyseMuon(ExRootTreeBranch *branch);
  void AnalyseTau(ExRootTreeBranch *branch);
  void AnalyseJet(ExRootTreeBranch *branch);
  void AnalyseMissingET(ExRootTreeBranch *branch);

  enum {kIntParamSize = 2, kDblParamSize = 7};
  Int_t fIntParam[kIntParamSize];
  Double_t fDblParam[kDblParamSize];

  Bool_t fIsReadyToFill;

  Int_t fTriggerWord, fEventNumber;

  char *fBuffer;

  ExRootTreeWriter *fTreeWriter;

  ExRootTreeBranch *fBranchEvent;
  ExRootTreeBranch *fBranchTrack;
  ExRootTreeBranch *fBranchTower;
  ExRootTreeBranch *fBranchPhoton;
  ExRootTreeBranch *fBranchElectron;
  ExRootTreeBranch *fBranchMuon;
  ExRootTreeBranch *fBranchJet;
  ExRootTreeBranch *fBranchMissingET;

};

//------------------------------------------------------------------------------

LHCOConverter::LHCOConverter(TFile *outputFile) :
  fIsReadyToFill(kFALSE),
  fTriggerWord(0), fEventNumber(1),
  fBuffer(0), fTreeWriter(0)
{
  fBuffer = new char[kBufferSize];
  fTreeWriter = new ExRootTreeWriter(outputFile, "Delphes");

  // information about reconstructed event
  fBranchEvent = fTreeWriter->NewBranch("Event", LHCOEvent::Class());
  // reconstructed tracks
  fBranchTrack = fTreeWriter->NewBranch("Track", Track::Class());
  /// calorimeter towers
  fBranchTower = fTreeWriter->NewBranch("Tower", Tower::Class());
  // reconstructed photons
  fBranchPhoton = fTreeWriter->NewBranch("Photon", Photon::Class());
  // reconstructed electrons
  fBranchElectron = fTreeWriter->NewBranch("Electron", Electron::Class());
  // reconstructed muons
  fBranchMuon = fTreeWriter->NewBranch("Muon", Muon::Class());
  // reconstructed jets
  fBranchJet = fTreeWriter->NewBranch("Jet", Jet::Class());
  // missing transverse energy
  fBranchMissingET = fTreeWriter->NewBranch("MissingET", MissingET::Class());
}

//------------------------------------------------------------------------------

LHCOConverter::~LHCOConverter()
{
  if(fTreeWriter) delete fTreeWriter;
  if(fBuffer) delete[] fBuffer;
}

//------------------------------------------------------------------------------

Bool_t LHCOConverter::ReadLine(FILE *inputFile)
{
  int rc;

  if(!fgets(fBuffer, kBufferSize, inputFile)) return kFALSE;

  DelphesStream bufferStream(fBuffer);

  rc = bufferStream.ReadInt(fIntParam[0]);

  if(!rc)
  {
    return kTRUE;
  }

  if(fIntParam[0] == 0)
  {
    rc = bufferStream.ReadInt(fEventNumber)
      && bufferStream.ReadInt(fTriggerWord);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid event format" << endl;
      return kFALSE;
    }

    if(fIsReadyToFill && fTreeWriter)
    {
      fTreeWriter->Fill();
      fTreeWriter->Clear();
    }

    AnalyseEvent(fBranchEvent);
    fIsReadyToFill = kTRUE;
  }
  else
  {
    rc = bufferStream.ReadInt(fIntParam[1])
      && bufferStream.ReadDbl(fDblParam[0])
      && bufferStream.ReadDbl(fDblParam[1])
      && bufferStream.ReadDbl(fDblParam[2])
      && bufferStream.ReadDbl(fDblParam[3])
      && bufferStream.ReadDbl(fDblParam[4])
      && bufferStream.ReadDbl(fDblParam[5])
      && bufferStream.ReadDbl(fDblParam[6]);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid object format" << endl;
      return kFALSE;
    }

    switch(fIntParam[1])
    {
      case 0: AnalysePhoton(fBranchPhoton); break;
      case 1: AnalyseElectron(fBranchElectron); break;
      case 2: AnalyseMuon(fBranchMuon); break;
      case 3: AnalyseTau(fBranchJet); break;
      case 4: AnalyseJet(fBranchJet); break;
      case 6: AnalyseMissingET(fBranchMissingET); break;
    }
  }

  return kTRUE;
}

//---------------------------------------------------------------------------

void LHCOConverter::Write()
{
  if(fIsReadyToFill && fTreeWriter) fTreeWriter->Fill();
  if(fTreeWriter) fTreeWriter->Write();
  fIsReadyToFill = kFALSE;
}

//---------------------------------------------------------------------------

void LHCOConverter::AnalyseEvent(ExRootTreeBranch *branch)
{
  LHCOEvent *element;

  element = static_cast<LHCOEvent*>(branch->NewEntry());

  element->Number = fEventNumber;
  element->Trigger = fTriggerWord;
}

//---------------------------------------------------------------------------

void LHCOConverter::AnalysePhoton(ExRootTreeBranch *branch)
{
  Photon *element;

  element = static_cast<Photon*>(branch->NewEntry());

  element->Eta = fDblParam[0];
  element->Phi = fDblParam[1];
  element->PT = fDblParam[2];
  element->EhadOverEem = fDblParam[6];
}

//---------------------------------------------------------------------------

void LHCOConverter::AnalyseElectron(ExRootTreeBranch *branch)
{
  Electron *element;

  element = static_cast<Electron*>(branch->NewEntry());

  element->Eta = fDblParam[0];
  element->Phi = fDblParam[1];
  element->PT = fDblParam[2];

  element->Charge = fDblParam[4] < 0.0 ? -1 : 1;
/*
  element->Ntrk = TMath::Abs(fDblParam[4]);
*/
  element->EhadOverEem = fDblParam[6];
}

//---------------------------------------------------------------------------

void LHCOConverter::AnalyseMuon(ExRootTreeBranch *branch)
{
  Muon *element;

  element = static_cast<Muon*>(branch->NewEntry());

  element->Eta = fDblParam[0];
  element->Phi = fDblParam[1];
  element->PT = fDblParam[2];

  element->Charge = fDblParam[4] < 0.0 ? -1 : 1;
/*
  element->Ntrk = TMath::Abs(fDblParam[4]);

  element->JetIndex = Int_t(fDblParam[5]);

  element->PTiso = Int_t(fDblParam[6]);
  element->ETiso = fDblParam[6] - element->PTiso;
*/
}

//---------------------------------------------------------------------------

void LHCOConverter::AnalyseTau(ExRootTreeBranch *branch)
{
  Jet *element;

  element = static_cast<Jet*>(branch->NewEntry());

  element->Eta = fDblParam[0];
  element->Phi = fDblParam[1];
  element->PT = fDblParam[2];

  element->Mass = fDblParam[3];

  element->BTag = 0;
  element->TauTag = 1;

  element->Charge = fDblParam[4] < 0 ? -1 : 1;
/*
  element->Ntrk = TMath::Abs(fDblParam[4]);
*/
  element->EhadOverEem = fDblParam[6];
}

//---------------------------------------------------------------------------

void LHCOConverter::AnalyseJet(ExRootTreeBranch *branch)
{
  Jet *element;

  element = static_cast<Jet*>(branch->NewEntry());

  element->Eta = fDblParam[0];
  element->Phi = fDblParam[1];
  element->PT = fDblParam[2];

  element->Mass = fDblParam[3];
/*
  element->Ntrk = TMath::Abs(Int_t(fDblParam[4]));
*/
  element->BTag = Int_t(fDblParam[5]);
  element->TauTag = 0;

  element->Charge = 0;

  element->EhadOverEem = fDblParam[6];
/*
  element->Index = fIntParam[0];
*/
}

//---------------------------------------------------------------------------

void LHCOConverter::AnalyseMissingET(ExRootTreeBranch *branch)
{
  MissingET *element;

  element = static_cast<MissingET*>(branch->NewEntry());

  element->Phi = fDblParam[1];
  element->MET = fDblParam[2];
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
  char appName[] = "lhco2root";
  stringstream message;
  FILE *inputFile = 0;
  TFile *outputFile = 0;
  LHCOConverter *converter = 0;
  Int_t i;
  Long64_t length, eventCounter;

  if(argc < 2)
  {
    cout << " Usage: " << appName << " output_file" << " [input_file(s)]" << endl;
    cout << " output_file - output file in ROOT format," << endl;
    cout << " input_file(s) - input file(s) in LHCO format," << endl;
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
    outputFile = TFile::Open(argv[1], "CREATE");

    if(outputFile == NULL)
    {
      message << "can't open " << argv[1];
      throw runtime_error(message.str());
    }

    converter = new LHCOConverter(outputFile);

    i = 2;
    do
    {
      if(interrupted) break;

      if(i == argc || strcmp(argv[i], "-") == 0)
      {
        cout << "** Reading standard input" << endl;
        inputFile = stdin;
        length = -1;
      }
      else
      {
        cout << "** Reading " << argv[i] << endl;
        inputFile = fopen(argv[i], "r");

        if(inputFile == NULL)
        {
          message << "can't open " << argv[i];
          throw runtime_error(message.str());
        }

        fseek(inputFile, 0L, SEEK_END);
        length = ftello(inputFile);
        fseek(inputFile, 0L, SEEK_SET);

        if(length <= 0)
        {
          fclose(inputFile);
          ++i;
          continue;
        }
      }

      eventCounter = 0;
      ExRootProgressBar progressBar(length);

      // Loop over all objects
      while(converter->ReadLine(inputFile) && !interrupted)
      {
        ++eventCounter;

        progressBar.Update(ftello(inputFile), eventCounter);
      }
      converter->Write();

      fseek(inputFile, 0L, SEEK_END);
      progressBar.Update(ftello(inputFile), eventCounter, kTRUE);
      progressBar.Finish();

      if(inputFile != stdin) fclose(inputFile);

      ++i;
    }
    while(i < argc);

    cout << "** Exiting..." << endl;

    delete converter;
    delete outputFile;

    return 0;
  }
  catch(runtime_error &e)
  {
    if(converter) delete converter;
    if(outputFile) delete outputFile;
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}


