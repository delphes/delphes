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

#include <stdexcept>
#include <iostream>
#include <sstream>

#include <stdio.h>
#include <errno.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "TObjArray.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "ExRootAnalysis/ExRootTreeBranch.h"

using namespace std;

static const int kBufferSize  = 1000000;

//---------------------------------------------------------------------------

DelphesSTDHEPReader::DelphesSTDHEPReader() :
  fInputFile(0), fInputXDR(0), fBuffer(0), fPDG(0), fBlockType(-1)
{
  fInputXDR = new XDR;
  fBuffer = new char[kBufferSize*96 + 24];

  fPDG = TDatabasePDG::Instance();
}

//---------------------------------------------------------------------------

DelphesSTDHEPReader::~DelphesSTDHEPReader()
{
  if(fBuffer) delete fBuffer;
  if(fInputXDR) delete fInputXDR;
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::SetInputFile(FILE *inputFile)
{
  fInputFile = inputFile;
  xdrstdio_create(fInputXDR, inputFile, XDR_DECODE);
  ReadFileHeader();
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
  if(feof(fInputFile)) return kFALSE;

  xdr_int(fInputXDR, &fBlockType);

  SkipBytes(4);

  if(fBlockType == EVENTTABLE)
  {
    ReadEventTable();
  }
  else if(fBlockType == EVENTHEADER)
  {
    ReadEventHeader();
  }
  else if(fBlockType == MCFIO_STDHEPBEG ||
          fBlockType == MCFIO_STDHEPEND)
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

void DelphesSTDHEPReader::SkipBytes(u_int size)
{
  int rc;
  u_int rndup;

  rndup = size % 4;
  if(rndup > 0)
  {
    rndup = 4 - rndup;
  }

  rc = fseek(fInputFile, size + rndup, SEEK_CUR);

  if(rc != 0 && errno == ESPIPE)
  {
    xdr_opaque(fInputXDR, fBuffer, size);
  }
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::SkipArray(u_int elsize)
{
  u_int size;
  xdr_u_int(fInputXDR, &size);
  SkipBytes(size*elsize);
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::ReadFileHeader()
{
  u_int i;
  enum STDHEPVersion {UNKNOWN, V1, V2, V21} version;

  xdr_int(fInputXDR, &fBlockType);
  if (fBlockType != FILEHEADER)
  {
    throw runtime_error("Header block not found. File is probably corrupted.");
  }

  SkipBytes(4);

  // version
  xdr_string(fInputXDR, &fBuffer, 100);
  if(fBuffer[0] == '\0' || fBuffer[1] == '\0') version = UNKNOWN;
  else if(fBuffer[0] == '1') version = V1;
  else if(strncmp(fBuffer, "2.01", 4) == 0) version = V21;
  else if(fBuffer[0] == '2') version = V2;
  else version = UNKNOWN;

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
  xdr_u_int(fInputXDR, &fEntries);

  SkipBytes(8);

  // Number of blocks
  u_int nBlocks = 0;
  xdr_u_int(fInputXDR, &nBlocks);

  // Number of NTuples
  u_int nNTuples = 0;
  if(version != V1)
  {
    xdr_u_int(fInputXDR, &nNTuples);
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
  xdr_string(fInputXDR, &fBuffer, 100);
  if(strncmp(fBuffer, "1.00", 4) == 0)
  {
    SkipBytes(8);

    SkipArray(4);
    SkipArray(4);
    SkipArray(4);
    SkipArray(4);
    SkipArray(4);
  }
  else if(strncmp(fBuffer, "2.00", 4) == 0)
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
  u_int skipSize = 4;

  // version
  xdr_string(fInputXDR, &fBuffer, 100);
  if(strncmp(fBuffer, "2.00", 4) == 0)
  {
    skipNTuples = true;
  }
  else if(strncmp(fBuffer, "3.00", 4) == 0)
  {
    skipNTuples = true;
    skipSize = 8;
  }

  SkipBytes(20);

  u_int dimBlocks = 0;
  xdr_u_int(fInputXDR, &dimBlocks);

  u_int dimNTuples = 0;
  if(skipNTuples)
  {
    SkipBytes(4);
    xdr_u_int(fInputXDR, &dimNTuples);
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
  xdr_string(fInputXDR, &fBuffer, 100);

  // skip 5*4 + 2*8 = 36 bytes
  SkipBytes(36);

  if((strncmp(fBuffer, "1.", 2) == 0) || (strncmp(fBuffer, "2.", 2) == 0) ||
     (strncmp(fBuffer, "3.", 2) == 0) || (strncmp(fBuffer, "4.", 2) == 0) ||
     (strncmp(fBuffer, "5.00", 4) == 0))
  {
    return;
  }

  SkipArray(1);
  SkipArray(1);

  if(strncmp(fBuffer, "5.01", 4) == 0)
  {
    return;
  }

  SkipBytes(4);
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::ReadSTDHEP()
{
  u_int idhepSize, isthepSize, jmohepSize, jdahepSize, phepSize, vhepSize;

  // version
  xdr_string(fInputXDR, &fBuffer, 100);

  // Extracting the event number
  xdr_int(fInputXDR, &fEventNumber);

  // Extracting the number of particles
  xdr_int(fInputXDR, &fEventSize);

  if(fEventSize >= kBufferSize)
  {
    throw runtime_error("too many particles in event");
  }

  // 4*n + 4*n + 8*n + 8*n + 40*n + 32*n +
  // 4 + 4 + 4 + 4 + 4 + 4 = 96*n + 24

  xdr_opaque(fInputXDR, fBuffer, 96*fEventSize + 24);

  idhepSize = ntohl(*(u_int*)(fBuffer));
  isthepSize = ntohl(*(u_int*)(fBuffer + 4*1 + 4*1*fEventSize));
  jmohepSize = ntohl(*(u_int*)(fBuffer + 4*2 + 4*2*fEventSize));
  jdahepSize = ntohl(*(u_int*)(fBuffer + 4*3 + 4*4*fEventSize));
  phepSize = ntohl(*(u_int*)(fBuffer + 4*4 + 4*6*fEventSize));
  vhepSize = ntohl(*(u_int*)(fBuffer + 4*5 + 4*16*fEventSize));

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
  memset(fScale, 0, 10*sizeof(double));
}

//---------------------------------------------------------------------------

void DelphesSTDHEPReader::ReadSTDHEP4()
{
  u_int number;

  // Extracting the event weight
  xdr_double(fInputXDR, &fWeight);

  // Extracting alpha QED
  xdr_double(fInputXDR, &fAlphaQED);

  // Extracting alpha QCD
  xdr_double(fInputXDR, &fAlphaQCD);

  // Extracting the event scale
  xdr_u_int(fInputXDR, &fScaleSize);
  for(number = 0; number < fScaleSize; ++number)
  {
    xdr_double(fInputXDR, &fScale[number]);
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
  int pid, status, m1, m2, d1, d2;
  double px, py, pz, e, mass;
  double x, y, z, t;

  XDR bufferXDR[6];
  xdrmem_create(&bufferXDR[0], fBuffer + 4*1, 4*fEventSize, XDR_DECODE);
  xdrmem_create(&bufferXDR[1], fBuffer + 4*2 + 4*1*fEventSize, 4*fEventSize, XDR_DECODE);
  xdrmem_create(&bufferXDR[2], fBuffer + 4*3 + 4*2*fEventSize, 8*fEventSize, XDR_DECODE);
  xdrmem_create(&bufferXDR[3], fBuffer + 4*4 + 4*4*fEventSize, 8*fEventSize, XDR_DECODE);
  xdrmem_create(&bufferXDR[4], fBuffer + 4*5 + 4*6*fEventSize, 40*fEventSize, XDR_DECODE);
  xdrmem_create(&bufferXDR[5], fBuffer + 4*6 + 4*16*fEventSize, 32*fEventSize, XDR_DECODE);

  for(number = 0; number < fEventSize; ++number)
  {
    xdr_int(&bufferXDR[0], &status);
    xdr_int(&bufferXDR[1], &pid);
    xdr_int(&bufferXDR[2], &m1);
    xdr_int(&bufferXDR[2], &m2);
    xdr_int(&bufferXDR[3], &d1);
    xdr_int(&bufferXDR[3], &d2);

    xdr_double(&bufferXDR[4], &px);
    xdr_double(&bufferXDR[4], &py);
    xdr_double(&bufferXDR[4], &pz);
    xdr_double(&bufferXDR[4], &e);
    xdr_double(&bufferXDR[4], &mass);

    xdr_double(&bufferXDR[5], &x);
    xdr_double(&bufferXDR[5], &y);
    xdr_double(&bufferXDR[5], &z);
    xdr_double(&bufferXDR[5], &t);

    candidate = factory->NewCandidate();

    candidate->PID = pid;
    pdgCode = TMath::Abs(candidate->PID);

    candidate->Status = status;

    candidate->M1 = m1 - 1;
    candidate->M2 = m2 - 1;

    candidate->D1 = d1 - 1;
    candidate->D2 = d2 - 1;

    pdgParticle = fPDG->GetParticle(pid);
    candidate->Charge = pdgParticle ? int(pdgParticle->Charge()/3.0) : -999;
    candidate->Mass = mass;

    candidate->Momentum.SetPxPyPzE(px, py, pz, e);

    candidate->Position.SetXYZT(x, y, z, t);

    allParticleOutputArray->Add(candidate);

    if(!pdgParticle) continue;

    if(status == 1 && pdgParticle->Stable())
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
