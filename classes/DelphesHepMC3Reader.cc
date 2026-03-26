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

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesHepMC3Reader.h"
#include "classes/DelphesStream.h"

#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TParticlePDG.h>
#include <TStopwatch.h>

#include <ExRootAnalysis/ExRootProgressBar.h>

#include <stdio.h>

DelphesHepMC3Reader::DelphesHepMC3Reader(const DelphesParameters &readerParams) :
  DelphesReader(readerParams),
  fPDG(TDatabasePDG::Instance()), fVertexCounter(-2), fParticleCounter(-1) {}

//---------------------------------------------------------------------------

void DelphesHepMC3Reader::LoadInputFile(std::string_view inputFile)
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

void DelphesHepMC3Reader::SetFactory(DelphesFactory *factory)
{
  DelphesReader::SetFactory(factory);
  fEventObject = GetFactory()->Book<HepMCEvent>("Event", true);
  fWeightsObject = GetFactory()->Book<std::vector<Weight> >("Weight", true);
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
  fWeightsObject->clear();
}

//---------------------------------------------------------------------------

bool DelphesHepMC3Reader::EventReady()
{
  return (fVertexCounter == -1) && (fParticleCounter == 0);
}

//---------------------------------------------------------------------------

bool DelphesHepMC3Reader::ReadBlock()
{
  std::map<int, std::pair<int, int> >::iterator itDaughterMap;
  char key, momentumUnit[4], positionUnit[3];
  int rc, code;
  double weight;

  if(!fgets(fBuffer.data(), fBuffer.size(), fInputFile)) return false;

  DelphesStream bufferStream(fBuffer.data() + 1);

  key = fBuffer.at(0);

  if(key == 'E')
  {
    Clear();

    int eventNumber;

    rc = bufferStream.ReadInt(eventNumber)
      && bufferStream.ReadInt(fVertexCounter)
      && bufferStream.ReadInt(fParticleCounter);

    fEventObject->Number = eventNumber; // int -> long

    if(!rc)
    {
      std::cerr << "** ERROR: invalid event format" << std::endl;
      return false;
    }
  }
  else if(key == 'U')
  {
    rc = sscanf(fBuffer.data() + 1, "%3s %2s", momentumUnit, positionUnit);

    if(rc != 2)
    {
      std::cerr << "** ERROR: invalid units format" << std::endl;
      return false;
    }

    if(strncmp(momentumUnit, "GEV", 3) == 0)
      fMomentumCoefficient = 1.0;
    else if(strncmp(momentumUnit, "MEV", 3) == 0)
      fMomentumCoefficient = 0.001;

    if(strncmp(positionUnit, "MM", 3) == 0)
      fPositionCoefficient = 1.0;
    else if(strncmp(positionUnit, "CM", 3) == 0)
      fPositionCoefficient = 10.0;
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
    rc = bufferStream.ReadInt(fEventObject->MPI);

    if(!rc)
    {
      std::cerr << "** ERROR: invalid MPI format" << std::endl;
      return false;
    }
  }
  else if(key == 'A' && bufferStream.FindStr("signal_process_id"))
  {
    rc = bufferStream.ReadInt(fEventObject->ProcessID);

    if(!rc)
    {
      std::cerr << "** ERROR: invalid process ID format" << std::endl;
      return false;
    }
  }
  else if(key == 'A' && bufferStream.FindStr("event_scale"))
  {
    double scale;
    rc = bufferStream.ReadDbl(scale);
    fEventObject->Scale = scale; // double -> float

    if(!rc)
    {
      std::cerr << "** ERROR: invalid event scale format" << std::endl;
      return false;
    }
  }
  else if(key == 'A' && bufferStream.FindStr("alphaQCD"))
  {
    double alphaQCD;
    rc = bufferStream.ReadDbl(alphaQCD);
    fEventObject->AlphaQCD = alphaQCD; // double -> float

    if(!rc)
    {
      std::cerr << "** ERROR: invalid alphaQCD format" << std::endl;
      return false;
    }
  }
  else if(key == 'A' && bufferStream.FindStr("alphaQED"))
  {
    double alphaQED;
    rc = bufferStream.ReadDbl(alphaQED);
    fEventObject->AlphaQED = alphaQED; // double -> float

    if(!rc)
    {
      std::cerr << "** ERROR: invalid alphaQED format" << std::endl;
      return false;
    }
  }
  else if(key == 'A' && bufferStream.FindStr("GenCrossSection"))
  {
    double crossSection, crossSectionError;
    rc = bufferStream.ReadDbl(crossSection)
      && bufferStream.ReadDbl(crossSectionError);

    // double -> float
    fEventObject->CrossSection = crossSection;
    fEventObject->CrossSectionError = crossSectionError;

    if(!rc)
    {
      std::cerr << "** ERROR: invalid cross section format" << std::endl;
      return false;
    }
  }
  else if(key == 'A' && bufferStream.FindStr("GenPdfInfo"))
  {
    double scalePDF, x1, x2, pdf1, pdf2;

    rc = bufferStream.ReadInt(fEventObject->ID1)
      && bufferStream.ReadInt(fEventObject->ID2)
      && bufferStream.ReadDbl(x1)
      && bufferStream.ReadDbl(x2)
      && bufferStream.ReadDbl(scalePDF)
      && bufferStream.ReadDbl(pdf1)
      && bufferStream.ReadDbl(pdf2);

    // double -> float
    fEventObject->X1 = x1;
    fEventObject->X2 = x2;
    fEventObject->ScalePDF = scalePDF;
    fEventObject->PDF1 = pdf1;
    fEventObject->PDF2 = pdf2;

    if(!rc)
    {
      std::cerr << "** ERROR: invalid PDF format" << std::endl;
      return false;
    }
  }
  else if(key == 'V')
  {
    int vertexCode, vertexStatus;
    fParticles.clear();

    fX = 0.0;
    fY = 0.0;
    fZ = 0.0;
    fT = 0.0;

    rc = bufferStream.ReadInt(vertexCode)
      && bufferStream.ReadInt(vertexStatus);

    if(!rc)
    {
      std::cerr << "** ERROR: invalid vertex format" << std::endl;
      return false;
    }

    rc = bufferStream.FindChr('[');

    if(!rc)
    {
      std::cerr << "** ERROR: invalid vertex format" << std::endl;
      return false;
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
        std::cerr << "** ERROR: invalid vertex format" << std::endl;
        return false;
      }
    }

    AnalyzeVertex(vertexCode);
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
      std::cerr << "** ERROR: invalid particle format" << std::endl;
      return false;
    }

    AnalyzeParticle();
  }

  if(EventReady())
  {
    FinalizeParticles();
  }

  return true;
}

//---------------------------------------------------------------------------

void DelphesHepMC3Reader::AnalyzeEvent(TStopwatch *procStopWatch)
{
  HepMCEvent &element = *fEventObject;

  element.Weight = fWeights.size() > 0 ? fWeights[0] : 1.0;

  element.ReadTime = fReadStopWatch.RealTime();
  element.ProcTime = procStopWatch->RealTime();

  for(std::vector<double>::const_iterator itWeights = fWeights.begin(); itWeights != fWeights.end(); ++itWeights)
  {
    Weight &weight = fWeightsObject->emplace_back();
    weight.Weight = *itWeights;
  }
}

//---------------------------------------------------------------------------

void DelphesHepMC3Reader::AnalyzeVertex(int code, Candidate *candidate)
{
  int index;
  std::shared_ptr<TLorentzVector> position;
  CandidatesCollection array;
  std::vector<int>::iterator itParticle;
  std::map<int, int>::iterator itVertexMap;

  itVertexMap = fOutVertexMap.find(code);
  if(itVertexMap == fOutVertexMap.end())
  {
    --fVertexCounter;

    index = fVertices.size();
    fOutVertexMap[code] = index;
    if(candidate && code > 0) fInVertexMap[code] = index;

    position = std::make_shared<TLorentzVector>();
    array = std::make_shared<std::vector<Candidate *> >();
    position->SetXYZT(0.0, 0.0, 0.0, 0.0);
    fVertices.push_back(std::make_pair(position, array));
  }
  else
  {
    index = itVertexMap->second;
    position = fVertices[index].first;
    array = fVertices[index].second;
  }

  if(candidate)
  {
    array->emplace_back(candidate);
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

void DelphesHepMC3Reader::AnalyzeParticle()
{
  Candidate *candidate = GetFactory()->NewCandidate();

  candidate->PID = fPID;
  candidate->Status = fParticleStatus;
  candidate->Mass = fMass;
  candidate->Momentum.SetPxPyPzE(fPx, fPy, fPz, fE);
  candidate->D1 = fParticleCode;

  AnalyzeVertex(fOutVertexCode, candidate);
}

//---------------------------------------------------------------------------

void DelphesHepMC3Reader::FinalizeParticles()
{
  TLorentzVector *position;
  CandidatesCollection array;
  Candidate *candidate;
  Candidate *candidateDaughter;
  TParticlePDG *pdgParticle;
  int pdgCode;
  std::map<int, int>::iterator itVertexMap;
  std::map<int, std::pair<int, int> >::iterator itMotherMap;
  std::map<int, std::pair<int, int> >::iterator itDaughterMap;
  int code, counter;

  counter = 0;
  for(size_t i = 0; i < fVertices.size(); ++i)
  {
    position = fVertices[i].first.get();
    array = fVertices[i].second;

    for(size_t j = 0; j < array->size(); ++j)
    {
      candidate = static_cast<Candidate *>(array->at(j));

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
        fDaughterMap[i] = std::make_pair(counter, counter);
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
          fMotherMap[code] = std::make_pair(counter, -1);
        }
        else
        {
          itMotherMap->second.second = counter;
        }
      }

      fAllParticleOutputArray->emplace_back(candidate);

      ++counter;

      pdgParticle = fPDG->GetParticle(candidate->PID);

      candidate->Charge = pdgParticle ? int(pdgParticle->Charge() / 3.0) : -999;

      if(!pdgParticle) continue;

      pdgCode = TMath::Abs(candidate->PID);

      if(candidate->Status == 1)
      {
        fStableParticleOutputArray->emplace_back(candidate);
      }
      else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15)
      {
        fPartonOutputArray->emplace_back(candidate);
      }
    }
  }

  for(size_t j = 0; j < fAllParticleOutputArray->size(); ++j)
  {
    candidate = static_cast<Candidate *>(fAllParticleOutputArray->at(j));

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

REGISTER_READER("HepMC3", DelphesHepMC3Reader);
