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


/** \class DelphesLHEFReader
 *
 *  Reads LHEF file
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
  fEventReady(kFALSE), fEventCounter(-1), fParticleCounter(-1)
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
  fEventReady = kFALSE;
  fEventCounter = -1;
  fParticleCounter = -1;
  fWeightList.clear();
}

//---------------------------------------------------------------------------

bool DelphesLHEFReader::EventReady()
{
  return fEventReady;
}

//---------------------------------------------------------------------------

bool DelphesLHEFReader::ReadBlock(DelphesFactory *factory,
  TObjArray *allParticleOutputArray,
  TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray)
{
  int rc, id;
  char *pch;
  double weight;

  if(!fgets(fBuffer, kBufferSize, fInputFile)) return kFALSE;

  if(strstr(fBuffer, "<event>"))
  {
    Clear();
    fEventCounter = 1;
  }
  else if(fEventCounter > 0)
  {
    DelphesStream bufferStream(fBuffer);

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
    DelphesStream bufferStream(fBuffer);

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
  else if(strstr(fBuffer, "<wgt"))
  {
    pch = strpbrk(fBuffer, "\"'");
    if(!pch)
    {
      cerr << "** ERROR: " << "invalid weight format" << endl;
      return kFALSE;
    }

    DelphesStream idStream(pch + 1);
    rc = idStream.ReadInt(id);

    pch = strchr(fBuffer, '>');
    if(!pch)
    {
      cerr << "** ERROR: " << "invalid weight format" << endl;
      return kFALSE;
    }

    DelphesStream weightStream(pch + 1);
    rc = weightStream.ReadDbl(weight);

    if(!rc)
    {
      cerr << "** ERROR: " << "invalid weight format" << endl;
      return kFALSE;
    }

    fWeightList.push_back(make_pair(id, weight));
  }
  else if(strstr(fBuffer, "</event>"))
  {
    fEventReady = kTRUE;
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

void DelphesLHEFReader::AnalyzeWeight(ExRootTreeBranch *branch)
{
  LHEFWeight *element;
  vector< pair< int, double > >::const_iterator itWeightList;

  for(itWeightList = fWeightList.begin(); itWeightList != fWeightList.end(); ++itWeightList)
  {
    element = static_cast<LHEFWeight *>(branch->NewEntry());

    element->ID = itWeightList->first;
    element->Weight = itWeightList->second;
  }
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
