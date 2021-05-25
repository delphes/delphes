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
  fVertexCounter(-1), fParticleCounter(-1)
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
  fWeight.clear();
  fMomentumCoefficient = 1.0;
  fPositionCoefficient = 1.0;
  fVertexCounter = -1;
  fParticleCounter = -1;
  fVertexMap.clear();
  fDaughterMap.clear();
}

//---------------------------------------------------------------------------

bool DelphesHepMC3Reader::EventReady()
{
  return (fParticleCounter == 0);
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

    fX = 0.0;
    fY = 0.0;
    fZ = 0.0;
    fT = 0.0;

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
      fWeight.push_back(weight);
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
    fX = 0.0;
    fY = 0.0;
    fZ = 0.0;
    fT = 0.0;

    fM1 = 0;
    fM2 = 0;

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
      if(code < fM1 || fM1 == 0) fM1 = code;
      if(code > fM2) fM2 = code;
      fVertexMap[code] = fVertexCode;
    }

    if(fM1 == fM2) fM2 = 0;

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

    itDaughterMap = fDaughterMap.find(fOutVertexCode);
    if(itDaughterMap == fDaughterMap.end())
    {
      fDaughterMap[fOutVertexCode] = make_pair(fParticleCode, fParticleCode);
    }
    else
    {
      itDaughterMap->second.second = fParticleCode;
    }

    AnalyzeParticle(factory, allParticleOutputArray,
      stableParticleOutputArray, partonOutputArray);
  }

  if(EventReady())
  {
    FinalizeParticles(allParticleOutputArray);
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
  element->Weight = fWeight.size() > 0 ? fWeight[0] : 1.0;
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
  vector<double>::const_iterator itWeight;

  for(itWeight = fWeight.begin(); itWeight != fWeight.end(); ++itWeight)
  {
    element = static_cast<Weight *>(branch->NewEntry());

    element->Weight = *itWeight;
  }
}

//---------------------------------------------------------------------------

void DelphesHepMC3Reader::AnalyzeParticle(DelphesFactory *factory,
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

  candidate->Status = fParticleStatus;

  pdgParticle = fPDG->GetParticle(fPID);
  candidate->Charge = pdgParticle ? int(pdgParticle->Charge() / 3.0) : -999;
  candidate->Mass = fMass;

  candidate->Momentum.SetPxPyPzE(fPx, fPy, fPz, fE);
  if(fMomentumCoefficient != 1.0)
  {
    candidate->Momentum *= fMomentumCoefficient;
  }

  candidate->Position.SetXYZT(fX, fY, fZ, fT);
  if(fPositionCoefficient != 1.0)
  {
    candidate->Position *= fPositionCoefficient;
  }

  candidate->D1 = -1;
  candidate->D2 = -1;

  if(fOutVertexCode < 0)
  {
    candidate->M1 = fM1 - 1;
    candidate->M2 = fM2 - 1;
  }
  else
  {
    candidate->M1 = fOutVertexCode - 1;
    candidate->M2 = -1;
  }

  allParticleOutputArray->Add(candidate);

  if(!pdgParticle) return;

  if(fParticleStatus == 1)
  {
    stableParticleOutputArray->Add(candidate);
  }
  else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15)
  {
    partonOutputArray->Add(candidate);
  }
}

//---------------------------------------------------------------------------

void DelphesHepMC3Reader::FinalizeParticles(TObjArray *allParticleOutputArray)
{
  Candidate *candidate;
  map<int, int >::iterator itVertexMap;
  map<int, pair<int, int> >::iterator itDaughterMap;
  int i, index;

  for(i = 0; i < allParticleOutputArray->GetEntriesFast(); ++i)
  {
    candidate = static_cast<Candidate *>(allParticleOutputArray->At(i));

    index = i + 1;

    itVertexMap = fVertexMap.find(index);
    if(itVertexMap != fVertexMap.end())
    {
      index = itVertexMap->second;
    }

    itDaughterMap = fDaughterMap.find(index);
    if(itDaughterMap == fDaughterMap.end())
    {
      candidate->D1 = -1;
      candidate->D2 = -1;
    }
    else
    {
      candidate->D1 = itDaughterMap->second.first - 1;
      candidate->D2 = itDaughterMap->second.second - 1;
    }
  }
}

//---------------------------------------------------------------------------
