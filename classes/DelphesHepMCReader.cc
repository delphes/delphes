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


/** \class DelphesHepMCReader
 *
 *  Reads HepMC file
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesHepMCReader.h"

#include <stdexcept>
#include <iostream>
#include <sstream>

#include <map>
#include <vector>

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
  fStateSize = 0;
  fState.clear();
  fWeightSize = 0;
  fWeight.clear();
  fMomentumCoefficient = 1.0;
  fPositionCoefficient = 1.0;
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
  char key, momentumUnit[4], positionUnit[3];
  int i, rc, state;
  double weight;

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
      && bufferStream.ReadInt(fVertexCounter)
      && bufferStream.ReadInt(fBeamCode[0])
      && bufferStream.ReadInt(fBeamCode[1])
      && bufferStream.ReadInt(fStateSize);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid event format" << endl;
      return kFALSE;
    }

    for(i = 0; i < fStateSize; ++i)
    {
      rc = rc && bufferStream.ReadInt(state);
      fState.push_back(state);
    }

    rc = rc && bufferStream.ReadInt(fWeightSize);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid event format" << endl;
      return kFALSE;
    }

    for(i = 0; i < fWeightSize; ++i)
    {
      rc = rc && bufferStream.ReadDbl(weight);
      fWeight.push_back(weight);
    }

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid event format" << endl;
      return kFALSE;
    }
  }
  else if(key == 'U')
  {
    rc = sscanf(fBuffer + 1, "%3s %2s", momentumUnit, positionUnit);

    if(rc != 2)
    {
      cerr << "** ERROR: " << "invalid units format" << endl;
      return kFALSE;
    }

    if(strncmp(momentumUnit, "GEV", 3) == 0)
    {
      fMomentumCoefficient = 1.0;
    }
    else if(strncmp(momentumUnit, "MEV", 3) == 0)
    {
      fMomentumCoefficient = 0.001;
    }
    
    if(strncmp(positionUnit, "MM", 3) == 0)
    {
      fPositionCoefficient = 1.0;
    }
    else if(strncmp(positionUnit, "CM", 3) == 0)
    {
      fPositionCoefficient = 10.0;
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
  element->Weight = fWeight.size() > 0 ? fWeight[0] : 1.0;
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
  candidate->Mass = fMass;

  candidate->Momentum.SetPxPyPzE(fPx, fPy, fPz, fE);
  if(fMomentumCoefficient != 1.0)
  {
    candidate->Momentum *= fMomentumCoefficient;
  }

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
    if(fPositionCoefficient != 1.0)
    {
      candidate->Position *= fPositionCoefficient;
    }
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
