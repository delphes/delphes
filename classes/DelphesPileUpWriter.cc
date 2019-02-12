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

/** \class DelphesPileUpWriter
 *
 *  Writes pile-up binary file
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesPileUpWriter.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <stdint.h>
#include <stdio.h>

#include "classes/DelphesXDRWriter.h"

using namespace std;

static const int kIndexSize = 10000000;
static const int kBufferSize = 1000000;
static const int kRecordSize = 9;

//------------------------------------------------------------------------------

DelphesPileUpWriter::DelphesPileUpWriter(const char *fileName) :
  fEntries(0), fEntrySize(0), fOffset(0),
  fPileUpFile(0), fIndex(0), fBuffer(0),
  fOutputWriter(0), fIndexWriter(0), fBufferWriter(0)
{
  stringstream message;

  fIndex = new uint8_t[kIndexSize * 8];
  fBuffer = new uint8_t[kBufferSize * kRecordSize * 4];
  fOutputWriter = new DelphesXDRWriter;
  fIndexWriter = new DelphesXDRWriter;
  fBufferWriter = new DelphesXDRWriter;

  fIndexWriter->SetBuffer(fIndex);
  fBufferWriter->SetBuffer(fBuffer);

  fPileUpFile = fopen(fileName, "wb");

  if(fPileUpFile == NULL)
  {
    message << "can't open pile-up file " << fileName;
    throw runtime_error(message.str());
  }

  fOutputWriter->SetFile(fPileUpFile);
}

//------------------------------------------------------------------------------

DelphesPileUpWriter::~DelphesPileUpWriter()
{
  if(fPileUpFile) fclose(fPileUpFile);
  if(fBufferWriter) delete fBufferWriter;
  if(fIndexWriter) delete fIndexWriter;
  if(fOutputWriter) delete fOutputWriter;
  if(fBuffer) delete[] fBuffer;
  if(fIndex) delete[] fIndex;
}

//------------------------------------------------------------------------------

void DelphesPileUpWriter::WriteParticle(int32_t pid,
  float x, float y, float z, float t,
  float px, float py, float pz, float e)
{
  if(fEntrySize >= kBufferSize)
  {
    throw runtime_error("too many particles in pile-up event");
  }

  fBufferWriter->WriteValue(&pid, 4);
  fBufferWriter->WriteValue(&x, 4);
  fBufferWriter->WriteValue(&y, 4);
  fBufferWriter->WriteValue(&z, 4);
  fBufferWriter->WriteValue(&t, 4);
  fBufferWriter->WriteValue(&px, 4);
  fBufferWriter->WriteValue(&py, 4);
  fBufferWriter->WriteValue(&pz, 4);
  fBufferWriter->WriteValue(&e, 4);

  ++fEntrySize;
}

//------------------------------------------------------------------------------

void DelphesPileUpWriter::WriteEntry()
{
  if(fEntries >= kIndexSize)
  {
    throw runtime_error("too many pile-up events");
  }

  fOutputWriter->WriteValue(&fEntrySize, 4);
  fOutputWriter->WriteRaw(fBuffer, fEntrySize * kRecordSize * 4);

  fIndexWriter->WriteValue(&fOffset, 8);
  fOffset += fEntrySize * kRecordSize * 4 + 4;

  fBufferWriter->SetOffset(0);
  fEntrySize = 0;

  ++fEntries;
}

//------------------------------------------------------------------------------

void DelphesPileUpWriter::WriteIndex()
{
  fOutputWriter->WriteRaw(fIndex, fEntries * 8);
  fOutputWriter->WriteValue(&fEntries, 8);
}

//------------------------------------------------------------------------------
