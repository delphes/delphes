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

/** \class DelphesSTDHEPReader
 *
 *  Reads STDHEP file
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesSTDHEPReader.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TParticlePDG.h"
#include "TStopwatch.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesXDRReader.h"

#include "ExRootAnalysis/ExRootTreeBranch.h"

using namespace std;

static const int kBufferSize = 1000000;

//---------------------------------------------------------------------------

DelphesSTDHEPReader::DelphesSTDHEPReader() :
  fInputFile(0), fBuffer(0), fPDG(0), fBlockType(-1)
{
  fBuffer = new uint8_t[kBufferSize * 96 + 24];

  fPDG = TDatabasePDG::Instance();
}

//---------------------------------------------------------------------------

DelphesSTDHEPReader::~DelphesSTDHEPReader()
{
  if(fBuffer) delete fBuffer;
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::SetInputFile(FILE *inputFile)
{
  fInputFile = inputFile;
  fReader[0].SetFile(inputFile);
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::Clear()
{
  fBlockType = -1;
}

//---------------------------------------------------------------------------

bool DelphesSTDHEPReader::EventReady()
{
  return (fBlockType == MCFIO_STDHEP) || (fBlockType == MCFIO_STDHEP4);
}

//---------------------------------------------------------------------------

bool DelphesSTDHEPReader::ReadBlock(DelphesFactory *factory,
  TObjArray *allParticleOutputArray,
  TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray)
{
  fReader[0].ReadValue(&fBlockType, 4);

  if(feof(fInputFile)) return kFALSE;

  SkipBytes(4);

  if(fBlockType == FILEHEADER)
  {
    ReadFileHeader();
  }
  else if(fBlockType == EVENTTABLE)
  {
    ReadEventTable();
  }
  else if(fBlockType == EVENTHEADER)
  {
    ReadEventHeader();
  }
  else if(fBlockType == MCFIO_STDHEPBEG || fBlockType == MCFIO_STDHEPEND)
  {
    ReadSTDCM1();
  }
  else if(fBlockType == MCFIO_STDHEP)
  {
    ReadSTDHEP();
    AnalyzeParticles(factory, allParticleOutputArray,
      stableParticleOutputArray, partonOutputArray);
  }
  else if(fBlockType == MCFIO_STDHEP4)
  {
    ReadSTDHEP();
    AnalyzeParticles(factory, allParticleOutputArray,
      stableParticleOutputArray, partonOutputArray);
    ReadSTDHEP4();
  }
  else
  {
    throw runtime_error("Unsupported block type.");
  }

  return kTRUE;
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::SkipBytes(int size)
{
  int rc;
  int rndup;

  rndup = size % 4;
  if(rndup > 0)
  {
    rndup = 4 - rndup;
  }

  rc = fseek(fInputFile, size + rndup, SEEK_CUR);

  if(rc != 0 && errno == ESPIPE)
  {
    fReader[0].ReadRaw(fBuffer, size);
  }
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::SkipArray(int elsize)
{
  uint32_t size;
  fReader[0].ReadValue(&size, 4);
  SkipBytes(size * elsize);
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::ReadFileHeader()
{
  uint32_t i;
  enum STDHEPVersion
  {
    UNKNOWN,
    V1,
    V2,
    V21
  } version;

  // version
  fReader[0].ReadString(fBuffer, 100);
  if(fBuffer[0] == '\0' || fBuffer[1] == '\0')
    version = UNKNOWN;
  else if(fBuffer[0] == '1')
    version = V1;
  else if(strncmp((char *)fBuffer, "2.01", 4) == 0)
    version = V21;
  else if(fBuffer[0] == '2')
    version = V2;
  else
    version = UNKNOWN;

  if(version == UNKNOWN)
  {
    throw runtime_error("Unknown file format version.");
  }

  SkipArray(1);
  SkipArray(1);
  SkipArray(1);

  if(version == V21)
  {
    SkipArray(1);
  }

  // Expected number of events
  SkipBytes(4);

  // Number of events
  fReader[0].ReadValue(&fEntries, 4);

  SkipBytes(8);

  // Number of blocks
  uint32_t nBlocks = 0;
  fReader[0].ReadValue(&nBlocks, 4);

  // Number of NTuples
  uint32_t nNTuples = 0;
  if(version != V1)
  {
    fReader[0].ReadValue(&nNTuples, 4);
  }

  if(nNTuples != 0)
  {
    throw runtime_error("Files containing n-tuples are not supported.");
  }

  // Processing blocks extraction
  if(nBlocks != 0)
  {
    SkipArray(4);

    for(i = 0; i < nBlocks; i++)
    {
      SkipArray(1);
    }
  }
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::ReadEventTable()
{
  // version
  fReader[0].ReadString(fBuffer, 100);
  if(strncmp((char *)fBuffer, "1.00", 4) == 0)
  {
    SkipBytes(8);

    SkipArray(4);
    SkipArray(4);
    SkipArray(4);
    SkipArray(4);
    SkipArray(4);
  }
  else if(strncmp((char *)fBuffer, "2.00", 4) == 0)
  {
    SkipBytes(12);

    SkipArray(4);
    SkipArray(4);
    SkipArray(4);
    SkipArray(4);
    SkipArray(8);
  }
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::ReadEventHeader()
{
  bool skipNTuples = false;
  int skipSize = 4;

  // version
  fReader[0].ReadString(fBuffer, 100);
  if(strncmp((char *)fBuffer, "2.00", 4) == 0)
  {
    skipNTuples = true;
  }
  else if(strncmp((char *)fBuffer, "3.00", 4) == 0)
  {
    skipNTuples = true;
    skipSize = 8;
  }

  SkipBytes(20);

  uint32_t dimBlocks = 0;
  fReader[0].ReadValue(&dimBlocks, 4);

  uint32_t dimNTuples = 0;
  if(skipNTuples)
  {
    SkipBytes(4);
    fReader[0].ReadValue(&dimNTuples, 4);
  }

  // Processing blocks extraction
  if(dimBlocks > 0)
  {
    SkipArray(4);
    SkipArray(skipSize);
  }

  // Processing blocks extraction
  if(skipNTuples && dimNTuples > 0)
  {
    SkipArray(4);
    SkipArray(skipSize);
  }
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::ReadSTDCM1()
{
  // version
  fReader[0].ReadString(fBuffer, 100);

  // skip 5*4 + 2*8 = 36 bytes
  SkipBytes(36);

  if((strncmp((char *)fBuffer, "1.", 2) == 0) || (strncmp((char *)fBuffer, "2.", 2) == 0) ||
     (strncmp((char *)fBuffer, "3.", 2) == 0) || (strncmp((char *)fBuffer, "4.", 2) == 0) ||
     (strncmp((char *)fBuffer, "5.00", 4) == 0))
  {
    return;
  }

  SkipArray(1);
  SkipArray(1);

  if(strncmp((char *)fBuffer, "5.01", 4) == 0)
  {
    return;
  }

  SkipBytes(4);
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::ReadSTDHEP()
{
  uint32_t idhepSize, isthepSize, jmohepSize, jdahepSize, phepSize, vhepSize;

  // version
  fReader[0].ReadString(fBuffer, 100);

  // Extracting the event number
  fReader[0].ReadValue(&fEventNumber, 4);

  // Extracting the number of particles
  fReader[0].ReadValue(&fEventSize, 4);

  if(fEventSize >= kBufferSize)
  {
    throw runtime_error("too many particles in event");
  }

  // 4*n + 4*n + 8*n + 8*n + 40*n + 32*n +
  // 4 + 4 + 4 + 4 + 4 + 4 = 96*n + 24

  fReader[0].ReadRaw(fBuffer, 96 * fEventSize + 24);

  fReader[1].SetBuffer(fBuffer);
  fReader[2].SetBuffer(fBuffer + 4 * 1 + 4 * 1 * fEventSize);
  fReader[3].SetBuffer(fBuffer + 4 * 2 + 4 * 2 * fEventSize);
  fReader[4].SetBuffer(fBuffer + 4 * 3 + 4 * 4 * fEventSize);
  fReader[5].SetBuffer(fBuffer + 4 * 4 + 4 * 6 * fEventSize);
  fReader[6].SetBuffer(fBuffer + 4 * 5 + 4 * 16 * fEventSize);

  fReader[1].ReadValue(&idhepSize, 4);
  fReader[2].ReadValue(&isthepSize, 4);
  fReader[3].ReadValue(&jmohepSize, 4);
  fReader[4].ReadValue(&jdahepSize, 4);
  fReader[5].ReadValue(&phepSize, 4);
  fReader[6].ReadValue(&vhepSize, 4);

  if(fEventSize < 0 ||
     fEventSize != (int)idhepSize      || fEventSize != (int)isthepSize     ||
     (2*fEventSize) != (int)jmohepSize || (2*fEventSize) != (int)jdahepSize ||
     (5*fEventSize) != (int)phepSize   || (4*fEventSize) != (int)vhepSize)
  {
    throw runtime_error("Inconsistent size of arrays. File is probably corrupted.");
  }

  fWeight = 1.0;
  fAlphaQED = 0.0;
  fAlphaQCD = 0.0;
  fScaleSize = 0;
  memset(fScale, 0, 10 * sizeof(double));
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::ReadSTDHEP4()
{
  uint32_t number;

  // Extracting the event weight
  fReader[0].ReadValue(&fWeight, 8);

  // Extracting alpha QED
  fReader[0].ReadValue(&fAlphaQED, 8);

  // Extracting alpha QCD
  fReader[0].ReadValue(&fAlphaQCD, 8);

  // Extracting the event scale
  fReader[0].ReadValue(&fScaleSize, 4);
  for(number = 0; number < fScaleSize; ++number)
  {
    fReader[0].ReadValue(&fScale[number], 8);
  }

  SkipArray(8);
  SkipArray(4);

  SkipBytes(4);
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::AnalyzeEvent(ExRootTreeBranch *branch, long long eventNumber,
  TStopwatch *readStopWatch, TStopwatch *procStopWatch)
{
  LHEFEvent *element;

  element = static_cast<LHEFEvent *>(branch->NewEntry());

  element->Number = fEventNumber;

  element->ProcessID = 0;

  element->Weight = fWeight;
  element->ScalePDF = fScale[0];
  element->AlphaQED = fAlphaQED;
  element->AlphaQCD = fAlphaQCD;

  element->ReadTime = readStopWatch->RealTime();
  element->ProcTime = procStopWatch->RealTime();
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::AnalyzeParticles(DelphesFactory *factory,
  TObjArray *allParticleOutputArray,
  TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray)
{
  Candidate *candidate;
  TParticlePDG *pdgParticle;
  int pdgCode;

  int number;
  int32_t pid, status, m1, m2, d1, d2;
  double px, py, pz, e, mass;
  double x, y, z, t;

  for(number = 0; number < fEventSize; ++number)
  {
    fReader[1].ReadValue(&status, 4);
    fReader[2].ReadValue(&pid, 4);
    fReader[3].ReadValue(&m1, 4);
    fReader[3].ReadValue(&m2, 4);
    fReader[4].ReadValue(&d1, 4);
    fReader[4].ReadValue(&d2, 4);

    fReader[5].ReadValue(&px, 8);
    fReader[5].ReadValue(&py, 8);
    fReader[5].ReadValue(&pz, 8);
    fReader[5].ReadValue(&e, 8);
    fReader[5].ReadValue(&mass, 8);

    fReader[6].ReadValue(&x, 8);
    fReader[6].ReadValue(&y, 8);
    fReader[6].ReadValue(&z, 8);
    fReader[6].ReadValue(&t, 8);

    candidate = factory->NewCandidate();

    candidate->PID = pid;
    pdgCode = TMath::Abs(candidate->PID);

    candidate->Status = status;

    candidate->M1 = m1 - 1;
    candidate->M2 = m2 - 1;

    candidate->D1 = d1 - 1;
    candidate->D2 = d2 - 1;

    pdgParticle = fPDG->GetParticle(pid);
    candidate->Charge = pdgParticle ? int(pdgParticle->Charge() / 3.0) : -999;
    candidate->Mass = mass;

    candidate->Momentum.SetPxPyPzE(px, py, pz, e);

    candidate->Position.SetXYZT(x, y, z, t);

    allParticleOutputArray->Add(candidate);

    if(!pdgParticle) continue;

    if(status == 1)
    {
      stableParticleOutputArray->Add(candidate);
    }
    else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15)
    {
      partonOutputArray->Add(candidate);
    }
  }
}

//---------------------------------------------------------------------------
