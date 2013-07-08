
/** \class DelphesLHEFReader
 *
 *  Reads LHEF file
 *
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesLHEFReader.h"

#include <stdexcept>
#include <iostream>
#include <sstream>

#include <stdio.h>

#include "TObjArray.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesStream.h"

#include "ExRootAnalysis/ExRootTreeBranch.h"

using namespace std;

static const int kBufferSize  = 1024;

//---------------------------------------------------------------------------

DelphesLHEFReader::DelphesLHEFReader() :
  fInputFile(0), fBuffer(0), fPDG(0),
  fEventCounter(-1), fParticleCounter(-1)
{
  fBuffer = new char[kBufferSize];

  fPDG = TDatabasePDG::Instance();
}

//---------------------------------------------------------------------------

DelphesLHEFReader::~DelphesLHEFReader()
{
  if(fBuffer) delete[] fBuffer;
}

//---------------------------------------------------------------------------

void DelphesLHEFReader::SetInputFile(FILE *inputFile)
{
  fInputFile = inputFile;
}

//---------------------------------------------------------------------------

void DelphesLHEFReader::Clear()
{
  fEventCounter = -1;
  fParticleCounter = -1;
}

//---------------------------------------------------------------------------

bool DelphesLHEFReader::EventReady()
{
  return (fEventCounter == 0) && (fParticleCounter == 0);
}

//---------------------------------------------------------------------------

bool DelphesLHEFReader::ReadBlock(DelphesFactory *factory,
  TObjArray *allParticleOutputArray,
  TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray)
{
  int rc;

  if(!fgets(fBuffer, kBufferSize, fInputFile)) return kFALSE;

  DelphesStream bufferStream(fBuffer);

  if(strstr(fBuffer, "<event>"))
  {
    Clear();
    fEventCounter = 1;
  }
  else if(fEventCounter > 0)
  {
    rc = bufferStream.ReadInt(fParticleCounter)
      && bufferStream.ReadInt(fProcessID)
      && bufferStream.ReadDbl(fWeight)
      && bufferStream.ReadDbl(fScalePDF)
      && bufferStream.ReadDbl(fAlphaQED)
      && bufferStream.ReadDbl(fAlphaQCD);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid event format" << endl;
      return kFALSE;
    }

    --fEventCounter;
  }
  else if(fParticleCounter > 0)
  {
    rc = bufferStream.ReadInt(fPID)
      && bufferStream.ReadInt(fStatus)
      && bufferStream.ReadInt(fM1)
      && bufferStream.ReadInt(fM2)
      && bufferStream.ReadInt(fC1)
      && bufferStream.ReadInt(fC2)
      && bufferStream.ReadDbl(fPx)
      && bufferStream.ReadDbl(fPy)
      && bufferStream.ReadDbl(fPz)
      && bufferStream.ReadDbl(fE)
      && bufferStream.ReadDbl(fMass);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid particle format" << endl;
      return kFALSE;
    }

    AnalyzeParticle(factory, allParticleOutputArray,
      stableParticleOutputArray, partonOutputArray);

    --fParticleCounter;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------

void DelphesLHEFReader::AnalyzeEvent(ExRootTreeBranch *branch, long long eventNumber,
  TStopwatch *readStopWatch, TStopwatch *procStopWatch)
{
  LHEFEvent *element;

  element = static_cast<LHEFEvent *>(branch->NewEntry());
  element->Number = eventNumber;

  element->ProcessID = fProcessID;
  element->Weight = fWeight;
  element->ScalePDF = fScalePDF;
  element->AlphaQED = fAlphaQED;
  element->AlphaQCD = fAlphaQCD;

  element->ReadTime = readStopWatch->RealTime();
  element->ProcTime = procStopWatch->RealTime();
}

//---------------------------------------------------------------------------

void DelphesLHEFReader::AnalyzeParticle(DelphesFactory *factory,
  TObjArray *allParticleOutputArray,
  TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray)
{
  Candidate *candidate;
  TParticlePDG *pdgParticle;
  int pdgCode;

  candidate = factory->NewCandidate();

  candidate->PID = fPID;
  pdgCode = TMath::Abs(candidate->PID);

  candidate->Status = fStatus;

  pdgParticle = fPDG->GetParticle(fPID);
  candidate->Charge = pdgParticle ? int(pdgParticle->Charge()/3.0) : -999;
  candidate->Mass = fMass;

  candidate->Momentum.SetPxPyPzE(fPx, fPy, fPz, fE);
  candidate->Position.SetXYZT(0.0, 0.0, 0.0, 0.0);

  candidate->M1 = fM1 - 1;
  candidate->M2 = fM2 - 1;

  candidate->D1 = -1;
  candidate->D2 = -1;

  allParticleOutputArray->Add(candidate);

  if(!pdgParticle) return;

  if(fStatus == 1 && pdgParticle->Stable())
  {
    stableParticleOutputArray->Add(candidate);
  }
  else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15)
  {
    partonOutputArray->Add(candidate);
  }
}

//---------------------------------------------------------------------------
