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

/** \class DelphesHepMC2Reader
 *
 *  Reads HepMC file
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesHepMC2Reader.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <map>
#include <vector>

#include <stdio.h>

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TParticlePDG.h"
#include "TStopwatch.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesStream.h"

#include <ExRootAnalysis/ExRootProgressBar.h>
#include <ExRootAnalysis/ExRootTreeBranch.h>

using namespace std;

DelphesHepMC2Reader::DelphesHepMC2Reader(const DelphesParameters &readerParams) :
  DelphesReader(readerParams),
  fInputFile(0), fPDG(TDatabasePDG::Instance()),
  fVertexCounter(-1), fInCounter(-1), fOutCounter(-1),
  fParticleCounter(0) {}

//---------------------------------------------------------------------------

void DelphesHepMC2Reader::LoadInputFile(std::string_view inputFile)
{
  if(fInputFile) fclose(fInputFile); // unload previous streams
  Clear(); // clear all buffers
  if(fInputFile = fopen(std::string{inputFile}.data(), "r"); fInputFile == nullptr)
  {
    std::ostringstream message;
    message << "can't open " << std::string{inputFile};
    throw std::runtime_error(message.str());
  }
  fseek(fInputFile, 0L, SEEK_END);
  int length = ftello(fInputFile);
  fProgressBar = std::make_unique<ExRootProgressBar>(length);
  fseek(fInputFile, 0L, SEEK_SET);

  if(length <= 0)
    fclose(fInputFile);
  fEventCounter = 0;
}

//---------------------------------------------------------------------------

void DelphesHepMC2Reader::SetFactory(DelphesFactory *factory)
{
  DelphesReader::SetFactory(factory);
  fEventObject = GetFactory()->Book<HepMCEvent>("Event");
  fWeights = GetFactory()->Book<std::vector<Weight> >("Weight");
}

//---------------------------------------------------------------------------

void DelphesHepMC2Reader::Clear()
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
  fWeights->clear();
}

//---------------------------------------------------------------------------

bool DelphesHepMC2Reader::EventReady()
{
  return (fVertexCounter == 0) && (fInCounter == 0) && (fOutCounter == 0);
}

//---------------------------------------------------------------------------

bool DelphesHepMC2Reader::ReadBlock()
{
  map<int, pair<int, int> >::iterator itMotherMap;
  map<int, pair<int, int> >::iterator itDaughterMap;
  char key, momentumUnit[4], positionUnit[3];
  int i, rc, state;
  double weight;

  if(!fgets(fBuffer.data(), fBuffer.size(), fInputFile)) return kFALSE;

  DelphesStream bufferStream(fBuffer.data() + 1);

  key = fBuffer.at(0);

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
      cerr << "** ERROR: "
           << "invalid event format" << endl;
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
      cerr << "** ERROR: "
           << "invalid event format" << endl;
      return kFALSE;
    }

    for(i = 0; i < fWeightSize; ++i)
    {
      rc = rc && bufferStream.ReadDbl(weight);
      fWeight.push_back(weight);
    }

    if(!rc)
    {
      cerr << "** ERROR: "
           << "invalid event format" << endl;
      return kFALSE;
    }
  }
  else if(key == 'U')
  {
    rc = sscanf(fBuffer.data() + 1, "%3s %2s", momentumUnit, positionUnit);

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
  else if(key == 'C')
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
      cerr << "** ERROR: "
           << "invalid PDF format" << endl;
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
      cerr << "** ERROR: "
           << "invalid vertex format" << endl;
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
      cerr << "** ERROR: "
           << "invalid particle format" << endl;
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

    AnalyzeParticle();

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
    FinalizeParticles();
  }

  return kTRUE;
}

//---------------------------------------------------------------------------

void DelphesHepMC2Reader::AnalyzeEvent(TStopwatch *procStopWatch)
{
  HepMCEvent &element = *fEventObject;

  element.Number = fEventNumber;

  element.ProcessID = fProcessID;
  element.MPI = fMPI;
  element.Weight = fWeight.size() > 0 ? fWeight[0] : 1.0;
  element.CrossSection = fCrossSection;
  element.CrossSectionError = fCrossSectionError;
  element.Scale = fScale;
  element.AlphaQED = fAlphaQED;
  element.AlphaQCD = fAlphaQCD;

  element.ID1 = fID1;
  element.ID2 = fID2;
  element.X1 = fX1;
  element.X2 = fX2;
  element.ScalePDF = fScalePDF;
  element.PDF1 = fPDF1;
  element.PDF2 = fPDF2;

  element.ReadTime = fReadStopWatch.RealTime();
  element.ProcTime = procStopWatch->RealTime();

  for(std::vector<double>::const_iterator itWeight = fWeight.begin(); itWeight != fWeight.end(); ++itWeight)
  {
    Weight &weight = fWeights->emplace_back();
    weight.Weight = *itWeight;
  }
}

//---------------------------------------------------------------------------

void DelphesHepMC2Reader::AnalyzeParticle()
{
  Candidate *candidate;
  TParticlePDG *pdgParticle;
  int pdgCode;

  candidate = GetFactory()->NewCandidate();

  candidate->PID = fPID;
  pdgCode = TMath::Abs(candidate->PID);

  candidate->Status = fStatus;

  pdgParticle = fPDG->GetParticle(fPID);
  candidate->Charge = pdgParticle ? int(pdgParticle->Charge() / 3.0) : -999;
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

  fAllParticleOutputArray->emplace_back(candidate);

  if(!pdgParticle) return;

  if(fStatus == 1)
  {
    fStableParticleOutputArray->emplace_back(candidate);
  }
  else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15)
  {
    fPartonOutputArray->emplace_back(candidate);
  }
}

//---------------------------------------------------------------------------

void DelphesHepMC2Reader::FinalizeParticles()
{
  Candidate *candidate;
  Candidate *candidateDaughter;
  map<int, pair<int, int> >::iterator itMotherMap;
  map<int, pair<int, int> >::iterator itDaughterMap;

  for(size_t i = 0; i < fAllParticleOutputArray->size(); ++i)
  {
    candidate = static_cast<Candidate *>(fAllParticleOutputArray->at(i));

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
        const TLorentzVector &decayPosition = candidate->Position;
        candidate->DecayPosition.SetXYZT(decayPosition.X(), decayPosition.Y(), decayPosition.Z(), decayPosition.T()); // decay position
      }
      else
      {
        candidate->D1 = itDaughterMap->second.first;
        candidate->D2 = itDaughterMap->second.second;
        candidateDaughter = static_cast<Candidate *>(fAllParticleOutputArray->at(candidate->D1));
        const TLorentzVector &decayPosition = candidateDaughter->Position;
        candidate->DecayPosition.SetXYZT(decayPosition.X(), decayPosition.Y(), decayPosition.Z(), decayPosition.T()); // decay position
      }
    }
  }
}

//---------------------------------------------------------------------------

REGISTER_READER("HepMC2", DelphesHepMC2Reader);
