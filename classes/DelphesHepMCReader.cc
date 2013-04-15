
/** \class DelphesHepMCReader
 *
 *  Reads HepMC file
 *
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesHepMCReader.h"

#include <stdexcept>
#include <iostream>
#include <sstream>

#include <map>

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

DelphesHepMCReader::DelphesHepMCReader() :
  fInputFile(0), fBuffer(0), fPDG(0),
  fVertexCounter(-1), fInCounter(-1), fOutCounter(-1),
  fParticleCounter(0)
{
  fBuffer = new char[kBufferSize];

  fPDG = TDatabasePDG::Instance();
}

//---------------------------------------------------------------------------

DelphesHepMCReader::~DelphesHepMCReader()
{
  if(fBuffer) delete[] fBuffer;
}

//---------------------------------------------------------------------------

void DelphesHepMCReader::SetInputFile(FILE *inputFile)
{
  fInputFile = inputFile;
}

//---------------------------------------------------------------------------

void DelphesHepMCReader::Clear()
{
  fVertexCounter = -1;
  fInCounter = -1;
  fOutCounter = -1;
  fMotherMap.clear();
  fDaughterMap.clear();
  fParticleCounter = 0;
}

//---------------------------------------------------------------------------

bool DelphesHepMCReader::EventReady()
{
  return (fVertexCounter == 0) && (fInCounter == 0) && (fOutCounter == 0);
}

//---------------------------------------------------------------------------

bool DelphesHepMCReader::ReadBlock(DelphesFactory *factory,
  TObjArray *allParticleOutputArray,
  TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray)
{
  map< int, pair< int, int > >::iterator itMotherMap;
  map< int, pair< int, int > >::iterator itDaughterMap;
  char key;
  int rc;

  if(!fgets(fBuffer, kBufferSize, fInputFile)) return kFALSE;

  DelphesStream bufferStream(fBuffer + 1);

  key = fBuffer[0];

  if(key == 'E')
  {
    Clear();

    rc = bufferStream.ReadInt(fEventNumber)
      && bufferStream.ReadInt(fMPI)
      && bufferStream.ReadDbl(fScale)
      && bufferStream.ReadDbl(fAlphaQCD)
      && bufferStream.ReadDbl(fAlphaQED)
      && bufferStream.ReadInt(fProcessID)
      && bufferStream.ReadInt(fSignalCode)
      && bufferStream.ReadInt(fVertexCounter);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid event format" << endl;
      return kFALSE;
    }
  }
  else if(key == 'F')
  {
    rc = bufferStream.ReadInt(fID1)
      && bufferStream.ReadInt(fID2)
      && bufferStream.ReadDbl(fX1)
      && bufferStream.ReadDbl(fX2)
      && bufferStream.ReadDbl(fScalePDF)
      && bufferStream.ReadDbl(fPDF1)
      && bufferStream.ReadDbl(fPDF2);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid PDF format" << endl;
      return kFALSE;
    }
  }
  else if(key == 'V' && fVertexCounter > 0)
  {
    rc = bufferStream.ReadInt(fOutVertexCode)
      && bufferStream.ReadInt(fVertexID)
      && bufferStream.ReadDbl(fX)
      && bufferStream.ReadDbl(fY)
      && bufferStream.ReadDbl(fZ)
      && bufferStream.ReadDbl(fT)
      && bufferStream.ReadInt(fInCounter)
      && bufferStream.ReadInt(fOutCounter);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid vertex format" << endl;
      return kFALSE;
    }
    --fVertexCounter;
  }
  else if(key == 'P' && fOutCounter > 0)
  {
    rc = bufferStream.ReadInt(fParticleCode)
      && bufferStream.ReadInt(fPID)
      && bufferStream.ReadDbl(fPx)
      && bufferStream.ReadDbl(fPy)
      && bufferStream.ReadDbl(fPz)
      && bufferStream.ReadDbl(fE)
      && bufferStream.ReadDbl(fMass)
      && bufferStream.ReadInt(fStatus)
      && bufferStream.ReadDbl(fTheta)
      && bufferStream.ReadDbl(fPhi)
      && bufferStream.ReadInt(fInVertexCode);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid particle format" << endl;
      return kFALSE;
    }

    if(fInVertexCode < 0)
    {
      itMotherMap = fMotherMap.find(fInVertexCode);
      if(itMotherMap == fMotherMap.end())
      {
        fMotherMap[fInVertexCode] = make_pair(fParticleCounter, -1);
      }
      else
      {
        itMotherMap->second.second = fParticleCounter;
      }
    }

    if(fInCounter <= 0)
    {
      itDaughterMap = fDaughterMap.find(fOutVertexCode);
      if(itDaughterMap == fDaughterMap.end())
      {
        fDaughterMap[fOutVertexCode] = make_pair(fParticleCounter, fParticleCounter);
      }
      else
      {
        itDaughterMap->second.second = fParticleCounter;
      }
    }

    AnalyzeParticle(factory, allParticleOutputArray,
      stableParticleOutputArray, partonOutputArray);

    if(fInCounter > 0)
    {
      --fInCounter;
    }
    else
    {
      --fOutCounter;
    }

    ++fParticleCounter;
  }

  if(EventReady())
  {
    FinalizeParticles(allParticleOutputArray);
  }

  return kTRUE;
}

//---------------------------------------------------------------------------

void DelphesHepMCReader::AnalyzeEvent(ExRootTreeBranch *branch, long long eventNumber,
  TStopwatch *readStopWatch, TStopwatch *procStopWatch)
{
  HepMCEvent *element;

  element = static_cast<HepMCEvent *>(branch->NewEntry());
  element->Number = fEventNumber;

  element->ProcessID = fProcessID;
  element->MPI = fMPI;
  element->Scale = fScale;
  element->AlphaQED = fAlphaQED;
  element->AlphaQCD = fAlphaQCD;

  element->ID1 = fID1;
  element->ID2 = fID2;
  element->X1 = fX1;
  element->X2 = fX2;
  element->ScalePDF = fScalePDF;
  element->PDF1 = fPDF1;
  element->PDF2 = fPDF2;

  element->ReadTime = readStopWatch->RealTime();
  element->ProcTime = procStopWatch->RealTime();
}

//---------------------------------------------------------------------------

void DelphesHepMCReader::AnalyzeParticle(DelphesFactory *factory,
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
  candidate->Mass = pdgParticle ? pdgParticle->Mass() : -999.9;

  candidate->Momentum.SetPxPyPzE(fPx, fPy, fPz, fE);

  candidate->M2 = 1;
  candidate->D2 = 1;
  if(fInCounter > 0)
  {
    candidate->M1 = 1;
    candidate->Position.SetXYZT(0.0, 0.0, 0.0, 0.0);
  }
  else
  {
    candidate->M1 = fOutVertexCode;
    candidate->Position.SetXYZT(fX, fY, fZ, fT);
  }
  if(fInVertexCode < 0)
  {
    candidate->D1 = fInVertexCode;
  }
  else
  {
    candidate->D1 = 1;
  }

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

void DelphesHepMCReader::FinalizeParticles(TObjArray *allParticleOutputArray)
{
  Candidate *candidate;
  map< int, pair< int, int > >::iterator itMotherMap;
  map< int, pair< int, int > >::iterator itDaughterMap;
  int i;

  for(i = 0; i < allParticleOutputArray->GetEntriesFast(); ++i)
  {
    candidate = static_cast<Candidate *>(allParticleOutputArray->At(i));

    if(candidate->M1 > 0)
    {
      candidate->M1 = -1;
      candidate->M2 = -1;
    }
    else
    {
      itMotherMap = fMotherMap.find(candidate->M1);
      if(itMotherMap == fMotherMap.end())
      {
        candidate->M1 = -1;
        candidate->M2 = -1;
      }
      else
      {
        candidate->M1 = itMotherMap->second.first;
        candidate->M2 = itMotherMap->second.second;
      }
    }
    if(candidate->D1 > 0)
    {
      candidate->D1 = -1;
      candidate->D2 = -1;
    }
    else
    {
      itDaughterMap = fDaughterMap.find(candidate->D1);
      if(itDaughterMap == fDaughterMap.end())
      {
        candidate->D1 = -1;
        candidate->D2 = -1;
      }
      else
      {
        candidate->D1 = itDaughterMap->second.first;
        candidate->D2 = itDaughterMap->second.second;
      }
    }
  }
}

//---------------------------------------------------------------------------
