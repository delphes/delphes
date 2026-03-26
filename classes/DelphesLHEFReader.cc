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

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesLHEFReader.h"
#include "classes/DelphesStream.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <stdio.h>

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TParticlePDG.h"
#include "TStopwatch.h"

#include <ExRootAnalysis/ExRootProgressBar.h>
#include <ExRootAnalysis/ExRootTreeBranch.h>

using namespace std;

static const int kBufferSize = 16384;

//---------------------------------------------------------------------------

DelphesLHEFReader::DelphesLHEFReader(const DelphesParameters &readerParams) :
  DelphesReader(readerParams),
  fPDG(TDatabasePDG::Instance()),
  fEventReady(false), fEventCounter(-1), fParticleCounter(-1), fCrossSection(1) {}

//---------------------------------------------------------------------------

void DelphesLHEFReader::LoadInputFile(std::string_view inputFile)
{
  if(fInputFile) fclose(fInputFile); // unload previous streams
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

void DelphesLHEFReader::SetFactory(DelphesFactory *factory)
{
  DelphesReader::SetFactory(factory);
  fEventInfo = GetFactory()->Book<LHEFEvent>("Event", true);
  fWeightInfo = GetFactory()->Book<std::vector<LHEFWeight> >("Weights", true);
}

//---------------------------------------------------------------------------

void DelphesLHEFReader::Clear()
{
  fEventReady = false;
  fEventCounter = -1;
  fParticleCounter = -1;
  fWeightList.clear();
  fWeightInfo->clear();
}

//---------------------------------------------------------------------------

bool DelphesLHEFReader::EventReady()
{
  return fEventReady;
}

//---------------------------------------------------------------------------

bool DelphesLHEFReader::ReadBlock()
{
  int rc, id;
  char *pch;
  double weight, xsec;

  if(!fgets(fBuffer.data(), fBuffer.size(), fInputFile)) return false;

  if(strstr(fBuffer.data(), "<event>"))
  {
    Clear();
    fEventCounter = 1;
  }
  else if(fEventCounter > 0)
  {
    DelphesStream bufferStream(fBuffer.data());

    rc = bufferStream.ReadInt(fParticleCounter)
      && bufferStream.ReadInt(fProcessID)
      && bufferStream.ReadDbl(fWeight)
      && bufferStream.ReadDbl(fScalePDF)
      && bufferStream.ReadDbl(fAlphaQED)
      && bufferStream.ReadDbl(fAlphaQCD);

    if(!rc)
    {
      cerr << "** ERROR: "
           << "invalid event format" << endl;
      return false;
    }

    --fEventCounter;
  }
  else if(fParticleCounter > 0)
  {
    DelphesStream bufferStream(fBuffer.data());

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
      cerr << "** ERROR: "
           << "invalid particle format" << endl;
      return false;
    }

    AnalyzeParticle();

    --fParticleCounter;
  }
  else if(strstr(fBuffer.data(), "<wgt"))
  {
    pch = strpbrk(fBuffer.data(), "\"'");
    if(!pch)
    {
      cerr << "** ERROR: "
           << "invalid weight format" << endl;
      return false;
    }

    DelphesStream idStream(pch + 1);
    rc = idStream.ReadInt(id);

    pch = strchr(fBuffer.data(), '>');
    if(!pch)
    {
      cerr << "** ERROR: "
           << "invalid weight format" << endl;
      return false;
    }

    DelphesStream weightStream(pch + 1);
    rc = weightStream.ReadDbl(weight);

    if(!rc)
    {
      cerr << "** ERROR: "
           << "invalid weight format" << endl;
      return false;
    }

    fWeightList.push_back(make_pair(id, weight));
  }
  else if(strstr(fBuffer.data(), "<xsecinfo"))
  {
    pch = strstr(fBuffer.data(), "totxsec");
    if(!pch)
    {
      cerr << "** ERROR: "
           << "invalid cross section format" << endl;
      return false;
    }

    pch = strpbrk(pch + 1, "\"'");
    if(!pch)
    {
      cerr << "** ERROR: "
           << "invalid cross section format" << endl;
      return false;
    }

    DelphesStream xsecStream(pch + 1);
    rc = xsecStream.ReadDbl(xsec);

    if(!rc)
    {
      cerr << "** ERROR: "
           << "invalid cross section format" << endl;
      return false;
    }

    fCrossSection = xsec;
  }
  else if(strstr(fBuffer.data(), "</event>"))
  {
    fEventReady = true;
  }
  return true;
}

//---------------------------------------------------------------------------

void DelphesLHEFReader::AnalyzeEvent(TStopwatch *procStopWatch)
{
  LHEFEvent &element = *fEventInfo;
  element.Number = fEventCounter;

  element.ProcessID = fProcessID;
  element.Weight = fWeight;
  element.CrossSection = fCrossSection;

  element.ScalePDF = fScalePDF;
  element.AlphaQED = fAlphaQED;
  element.AlphaQCD = fAlphaQCD;

  element.ReadTime = fReadStopWatch.RealTime();
  element.ProcTime = procStopWatch->RealTime();

  for(std::vector<pair<int, double> >::const_iterator itWeightList = fWeightList.begin();
    itWeightList != fWeightList.end(); ++itWeightList)
  {
    LHEFWeight &weight = fWeightInfo->emplace_back();
    weight.ID = itWeightList->first;
    weight.Weight = itWeightList->second;
  }
}

//---------------------------------------------------------------------------

void DelphesLHEFReader::AnalyzeParticle()
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
  candidate->Position.SetXYZT(0.0, 0.0, 0.0, 0.0);

  candidate->M1 = fM1 - 1;
  candidate->M2 = fM2 - 1;

  candidate->D1 = -1;
  candidate->D2 = -1;

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

REGISTER_READER("LHEF", DelphesLHEFReader);
