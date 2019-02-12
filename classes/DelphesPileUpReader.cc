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

/** \class DelphesPileUpReader
 *
 *  Reads pile-up binary file
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesPileUpReader.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <stdint.h>
#include <stdio.h>

#include "classes/DelphesXDRReader.h"

using namespace std;

static const int kIndexSize = 10000000;
static const int kBufferSize = 1000000;
static const int kRecordSize = 9;

//------------------------------------------------------------------------------

DelphesPileUpReader::DelphesPileUpReader(const char *fileName) :
  fEntries(0), fEntrySize(0), fCounter(0),
  fPileUpFile(0), fIndex(0), fBuffer(0),
  fInputReader(0), fIndexReader(0), fBufferReader(0)
{
  stringstream message;

  fIndex = new uint8_t[kIndexSize * 8];
  fBuffer = new uint8_t[kBufferSize * kRecordSize * 4];
  fInputReader = new DelphesXDRReader;
  fIndexReader = new DelphesXDRReader;
  fBufferReader = new DelphesXDRReader;

  fIndexReader->SetBuffer(fIndex);
  fBufferReader->SetBuffer(fBuffer);

  fPileUpFile = fopen(fileName, "rb");

  if(fPileUpFile == NULL)
  {
    message << "can't open pile-up file " << fileName;
    throw runtime_error(message.str());
  }

  fInputReader->SetFile(fPileUpFile);

  // read number of events
  fseeko(fPileUpFile, -8, SEEK_END);
  fInputReader->ReadValue(&fEntries, 8);

  if(fEntries >= kIndexSize)
  {
    message << "too many events in pile-up file " << fileName;
    throw runtime_error(message.str());
  }

  // read index of events
  fseeko(fPileUpFile, -8 - 8 * fEntries, SEEK_END);
  fInputReader->ReadRaw(fIndex, fEntries * 8);
}

//------------------------------------------------------------------------------

DelphesPileUpReader::~DelphesPileUpReader()
{
  if(fPileUpFile) fclose(fPileUpFile);
  if(fBufferReader) delete fBufferReader;
  if(fIndexReader) delete fIndexReader;
  if(fInputReader) delete fInputReader;
  if(fBuffer) delete[] fBuffer;
  if(fIndex) delete[] fIndex;
}

//------------------------------------------------------------------------------

bool DelphesPileUpReader::ReadParticle(int32_t &pid,
  float &x, float &y, float &z, float &t,
  float &px, float &py, float &pz, float &e)
{
  if(fCounter >= fEntrySize) return false;

  fBufferReader->ReadValue(&pid, 4);
  fBufferReader->ReadValue(&x, 4);
  fBufferReader->ReadValue(&y, 4);
  fBufferReader->ReadValue(&z, 4);
  fBufferReader->ReadValue(&t, 4);
  fBufferReader->ReadValue(&px, 4);
  fBufferReader->ReadValue(&py, 4);
  fBufferReader->ReadValue(&pz, 4);
  fBufferReader->ReadValue(&e, 4);

  ++fCounter;

  return true;
}

//------------------------------------------------------------------------------

bool DelphesPileUpReader::ReadEntry(int64_t entry)
{
  int64_t offset;

  if(entry >= fEntries) return false;

  // read event position
  fIndexReader->SetOffset(8 * entry);
  fIndexReader->ReadValue(&offset, 8);

  // read event
  fseeko(fPileUpFile, offset, SEEK_SET);
  fInputReader->ReadValue(&fEntrySize, 4);

  if(fEntrySize >= kBufferSize)
  {
    throw runtime_error("too many particles in pile-up event");
  }

  fInputReader->ReadRaw(fBuffer, fEntrySize * kRecordSize * 4);
  fBufferReader->SetOffset(0);
  fCounter = 0;

  return true;
}

//------------------------------------------------------------------------------
