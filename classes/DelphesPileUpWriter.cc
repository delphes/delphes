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

//------------------------------------------------------------------------------

DelphesPileUpWriter::DelphesPileUpWriter(std::string_view fileName) :
  fOutputWriter(std::make_unique<DelphesXDRWriter>()),
  fIndexWriter(std::make_unique<DelphesXDRWriter>()),
  fBufferWriter(std::make_unique<DelphesXDRWriter>())
{
  fIndexWriter->SetBuffer(fIndex.data());
  fBufferWriter->SetBuffer(fBuffer.data());

  if(fPileUpFile = fopen(fileName.data(), "wb"); fPileUpFile == nullptr)
  {
    std::ostringstream message;
    message << "can't open pile-up file " << fileName;
    throw std::runtime_error(message.str());
  }
  fOutputWriter->SetFile(fPileUpFile);
}

//------------------------------------------------------------------------------

DelphesPileUpWriter::~DelphesPileUpWriter()
{
  if(fPileUpFile) fclose(fPileUpFile);
}

//------------------------------------------------------------------------------

void DelphesPileUpWriter::WriteParticle(int32_t pid,
  float x, float y, float z, float t,
  float px, float py, float pz, float e)
{
  if(static_cast<size_t>(fEntrySize) >= kBufferSize)
    throw std::runtime_error("too many particles in pile-up event");

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
  if(static_cast<size_t>(fEntries) >= kIndexSize)
    throw std::runtime_error("too many pile-up events");

  fOutputWriter->WriteValue(&fEntrySize, 4);
  fOutputWriter->WriteRaw(fBuffer.data(), fEntrySize * kRecordSize * 4);

  fIndexWriter->WriteValue(&fOffset, 8);
  fOffset += fEntrySize * kRecordSize * 4 + 4;

  fBufferWriter->SetOffset(0);
  fEntrySize = 0;

  ++fEntries;
}

//------------------------------------------------------------------------------

void DelphesPileUpWriter::WriteIndex()
{
  fOutputWriter->WriteRaw(fIndex.data(), fEntries * 8);
  fOutputWriter->WriteValue(&fEntries, 8);
}

//------------------------------------------------------------------------------
