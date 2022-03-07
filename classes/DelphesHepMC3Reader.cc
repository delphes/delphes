/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2021  Universite catholique de Louvain (UCL), Belgium
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

/** \class DelphesHepMC3Reader
 *
 *  Reads HepMC file
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesHepMC3Reader.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <map>
#include <vector>

#include <stdio.h>

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TParticlePDG.h"
#include "TStopwatch.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesStream.h"

#include "ExRootAnalysis/ExRootTreeBranch.h"

using namespace std;

static const int kBufferSize = 16384;

//---------------------------------------------------------------------------

DelphesHepMC3Reader::DelphesHepMC3Reader() :
  fInputFile(0), fBuffer(0), fPDG(0),
  fVertexCounter(-2), fParticleCounter(-1)
{
  fBuffer = new char[kBufferSize];

  fPDG = TDatabasePDG::Instance();
}

//---------------------------------------------------------------------------

DelphesHepMC3Reader::~DelphesHepMC3Reader()
{
  if(fBuffer) delete[] fBuffer;
}

//---------------------------------------------------------------------------

void DelphesHepMC3Reader::SetInputFile(FILE *inputFile)
{
  fInputFile = inputFile;
}

//---------------------------------------------------------------------------

void DelphesHepMC3Reader::Clear()
{
  fWeights.clear();
  fMomentumCoefficient = 1.0;
  fPositionCoefficient = 1.0;
  fVertexCounter = -2;
  fParticleCounter = -1;
  fVertices.clear();
  fParticles.clear();
  fInVertexMap.clear();
  fOutVertexMap.clear();
  fMotherMap.clear();
  fDaughterMap.clear();
}

//---------------------------------------------------------------------------

bool DelphesHepMC3Reader::EventReady()
{
  return (fVertexCounter == -1) && (fParticleCounter == 0);
}

//---------------------------------------------------------------------------

bool DelphesHepMC3Reader::ReadBlock(DelphesFactory *factory,
  TObjArray *allParticleOutputArray,
  TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray)
{
  map<int, pair<int, int> >::iterator itDaughterMap;
  char key, momentumUnit[4], positionUnit[3];
  int rc, code;
  double weight;

  if(!fgets(fBuffer, kBufferSize, fInputFile)) return kFALSE;

  DelphesStream bufferStream(fBuffer + 1);

  key = fBuffer[0];

  if(key == 'E')
  {
    Clear();

    rc = bufferStream.ReadInt(fEventNumber)
      && bufferStream.ReadInt(fVertexCounter)
      && bufferStream.ReadInt(fParticleCounter);

    if(!rc)
    {
      cerr << "** ERROR: "
           << "invalid event format" << endl;
      return kFALSE;
    }
  }
  else if(key == 'U')
  {
    rc = sscanf(fBuffer + 1, "%3s %2s", momentumUnit, positionUnit);

    if(rc != 2)
    {
      cerr << "** ERROR: "
           << "invalid units format" << endl;
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
  else if(key == 'W')
  {
    while(bufferStream.ReadDbl(weight))
    {
      fWeights.push_back(weight);
    }
  }
  else if(key == 'A' && bufferStream.FindStr("mpi"))
  {
    rc = bufferStream.ReadInt(fMPI);

    if(!rc)
    {
      cerr << "** ERROR: "
           << "invalid MPI format" << endl;
      return kFALSE;
    }
  }
  else if(key == 'A' && bufferStream.FindStr("signal_process_id"))
  {
    rc = bufferStream.ReadInt(fProcessID);

    if(!rc)
    {
      cerr << "** ERROR: "
           << "invalid process ID format" << endl;
      return kFALSE;
    }
  }
  else if(key == 'A' && bufferStream.FindStr("event_scale"))
  {
    rc = bufferStream.ReadDbl(fScale);

    if(!rc)
    {
      cerr << "** ERROR: "
           << "invalid event scale format" << endl;
      return kFALSE;
    }
  }
  else if(key == 'A' && bufferStream.FindStr("alphaQCD"))
  {
    rc = bufferStream.ReadDbl(fAlphaQCD);

    if(!rc)
    {
      cerr << "** ERROR: "
           << "invalid alphaQCD format" << endl;
      return kFALSE;
    }
  }
  else if(key == 'A' && bufferStream.FindStr("alphaQED"))
  {
    rc = bufferStream.ReadDbl(fAlphaQED);

    if(!rc)
    {
      cerr << "** ERROR: "
           << "invalid alphaQED format" << endl;
      return kFALSE;
    }
  }
  else if(key == 'A' && bufferStream.FindStr("GenCrossSection"))
  {
    rc = bufferStream.ReadDbl(fCrossSection)
      && bufferStream.ReadDbl(fCrossSectionError);

    if(!rc)
    {
      cerr << "** ERROR: "
           << "invalid cross section format" << endl;
      return kFALSE;
    }
  }
  else if(key == 'A' && bufferStream.FindStr("GenPdfInfo"))
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
      cerr << "** ERROR: "
           << "invalid PDF format" << endl;
      return kFALSE;
    }
  }
  else if(key == 'V')
  {
    fParticles.clear();

    fX = 0.0;
    fY = 0.0;
    fZ = 0.0;
    fT = 0.0;

    rc = bufferStream.ReadInt(fVertexCode)
      && bufferStream.ReadInt(fVertexStatus);

    if(!rc)
    {
      cerr << "** ERROR: "
           << "invalid vertex format" << endl;
      return kFALSE;
    }

    rc = bufferStream.FindChr('[');

    if(!rc)
    {
      cerr << "** ERROR: "
           << "invalid vertex format" << endl;
      return kFALSE;
    }

    while(bufferStream.ReadInt(code))
    {
      fParticles.push_back(code);
      bufferStream.FindChr(',');
    }

    if(bufferStream.FindChr('@'))
    {
      rc = bufferStream.ReadDbl(fX)
        && bufferStream.ReadDbl(fY)
        && bufferStream.ReadDbl(fZ)
        && bufferStream.ReadDbl(fT);

      if(!rc)
      {
        cerr << "** ERROR: "
             << "invalid vertex format" << endl;
        return kFALSE;
      }
    }

    AnalyzeVertex(factory, fVertexCode);
  }
  else if(key == 'P' && fParticleCounter > 0)
  {
    --fParticleCounter;

    rc = bufferStream.ReadInt(fParticleCode)
      && bufferStream.ReadInt(fOutVertexCode)
      && bufferStream.ReadInt(fPID)
      && bufferStream.ReadDbl(fPx)
      && bufferStream.ReadDbl(fPy)
      && bufferStream.ReadDbl(fPz)
      && bufferStream.ReadDbl(fE)
      && bufferStream.ReadDbl(fMass)
      && bufferStream.ReadInt(fParticleStatus);

    if(!rc)
    {
      cerr << "** ERROR: "
           << "invalid particle format" << endl;
      return kFALSE;
    }

    AnalyzeParticle(factory);
  }

  if(EventReady())
  {
    FinalizeParticles(allParticleOutputArray, stableParticleOutputArray, partonOutputArray);
  }

  return kTRUE;
}

//---------------------------------------------------------------------------

void DelphesHepMC3Reader::AnalyzeEvent(ExRootTreeBranch *branch, long long eventNumber,
  TStopwatch *readStopWatch, TStopwatch *procStopWatch)
{
  HepMCEvent *element;

  element = static_cast<HepMCEvent *>(branch->NewEntry());
  element->Number = fEventNumber;

  element->ProcessID = fProcessID;
  element->MPI = fMPI;
  element->Weight = fWeights.size() > 0 ? fWeights[0] : 1.0;
  element->CrossSection = fCrossSection;
  element->CrossSectionError = fCrossSectionError;
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

void DelphesHepMC3Reader::AnalyzeWeight(ExRootTreeBranch *branch)
{
  Weight *element;
  vector<double>::const_iterator itWeights;

  for(itWeights = fWeights.begin(); itWeights != fWeights.end(); ++itWeights)
  {
    element = static_cast<Weight *>(branch->NewEntry());

    element->Weight = *itWeights;
  }
}

//---------------------------------------------------------------------------

void DelphesHepMC3Reader::AnalyzeVertex(DelphesFactory *factory, int code, Candidate *candidate)
{
  int index;
  TLorentzVector *position;
  TObjArray *array;
  vector<int>::iterator itParticle;
  map<int, int>::iterator itVertexMap;

  itVertexMap = fOutVertexMap.find(code);
  if(itVertexMap == fOutVertexMap.end())
  {
    --fVertexCounter;

    index = fVertices.size();
    fOutVertexMap[code] = index;
    if(candidate && code > 0) fInVertexMap[code] = index;

    position = factory->New<TLorentzVector>();
    array = factory->NewArray();
    position->SetXYZT(0.0, 0.0, 0.0, 0.0);
    fVertices.push_back(make_pair(position, array));
  }
  else
  {
    index = itVertexMap->second;
    position = fVertices[index].first;
    array = fVertices[index].second;
  }

  if(candidate)
  {
    array->Add(candidate);
  }
  else
  {
    position->SetXYZT(fX, fY, fZ, fT);
    for(itParticle = fParticles.begin(); itParticle != fParticles.end(); ++itParticle)
    {
      fInVertexMap[*itParticle] = index;
    }
  }
}

//---------------------------------------------------------------------------

void DelphesHepMC3Reader::AnalyzeParticle(DelphesFactory *factory)
{
  Candidate *candidate;

  candidate = factory->NewCandidate();

  candidate->PID = fPID;

  candidate->Status = fParticleStatus;

  candidate->Mass = fMass;

  candidate->Momentum.SetPxPyPzE(fPx, fPy, fPz, fE);

  candidate->D1 = fParticleCode;

  AnalyzeVertex(factory, fOutVertexCode, candidate);
}

//---------------------------------------------------------------------------

void DelphesHepMC3Reader::FinalizeParticles(TObjArray *allParticleOutputArray,
  TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray)
{
  TLorentzVector *position;
  TObjArray *array;
  Candidate *candidate;
  Candidate *candidateDaughter;
  TParticlePDG *pdgParticle;
  int pdgCode;
  map<int, int >::iterator itVertexMap;
  map<int, pair<int, int> >::iterator itMotherMap;
  map<int, pair<int, int> >::iterator itDaughterMap;
  int i, j, code, counter;

  counter = 0;
  for(i = 0; i < fVertices.size(); ++i)
  {
    position = fVertices[i].first;
    array = fVertices[i].second;

    for(j = 0; j < array->GetEntriesFast(); ++j)
    {
      candidate = static_cast<Candidate *>(array->At(j));

      candidate->Position = *position;
      if(fPositionCoefficient != 1.0)
      {
        candidate->Position *= fPositionCoefficient;
      }

      if(fMomentumCoefficient != 1.0)
      {
        candidate->Momentum *= fMomentumCoefficient;
      }

      candidate->M1 = i;

      itDaughterMap = fDaughterMap.find(i);
      if(itDaughterMap == fDaughterMap.end())
      {
        fDaughterMap[i] = make_pair(counter, counter);
      }
      else
      {
        itDaughterMap->second.second = counter;
      }

      code = candidate->D1;

      itVertexMap = fInVertexMap.find(code);
      if(itVertexMap == fInVertexMap.end())
      {
        candidate->D1 = -1;
      }
      else
      {
        code = itVertexMap->second;

        candidate->D1 = code;

        itMotherMap = fMotherMap.find(code);
        if(itMotherMap == fMotherMap.end())
        {
          fMotherMap[code] = make_pair(counter, -1);
        }
        else
        {
          itMotherMap->second.second = counter;
        }
      }

      allParticleOutputArray->Add(candidate);

      ++counter;

      pdgParticle = fPDG->GetParticle(candidate->PID);

      candidate->Charge = pdgParticle ? int(pdgParticle->Charge() / 3.0) : -999;

      if(!pdgParticle) continue;

      pdgCode = TMath::Abs(candidate->PID);

      if(candidate->Status == 1)
      {
        stableParticleOutputArray->Add(candidate);
      }
      else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15)
      {
        partonOutputArray->Add(candidate);
      }
    }
  }

  for(i = 0; i < allParticleOutputArray->GetEntriesFast(); ++i)
  {
    candidate = static_cast<Candidate *>(allParticleOutputArray->At(i));

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

    if(candidate->D1 < 0)
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
        const TLorentzVector &decayPosition = candidate->Position;
        candidate->DecayPosition.SetXYZT(decayPosition.X(), decayPosition.Y(), decayPosition.Z(), decayPosition.T());// decay position
      }
      else
      {
        candidate->D1 = itDaughterMap->second.first;
        candidate->D2 = itDaughterMap->second.second;
        candidateDaughter = static_cast<Candidate *>(allParticleOutputArray->At(candidate->D1));
        const TLorentzVector &decayPosition = candidateDaughter->Position;
        candidate->DecayPosition.SetXYZT(decayPosition.X(), decayPosition.Y(), decayPosition.Z(), decayPosition.T());// decay position
      }
    }
  }
}

//---------------------------------------------------------------------------
