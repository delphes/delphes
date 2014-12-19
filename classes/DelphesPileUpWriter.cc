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

#include <stdexcept>
#include <iostream>
#include <sstream>

#include <stdio.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

using namespace std;

static const int kIndexSize = 10000000;
static const int kBufferSize = 1000000;
static const int kRecordSize = 9;

//------------------------------------------------------------------------------

DelphesPileUpWriter::DelphesPileUpWriter(const char *fileName) :
  fEntries(0), fEntrySize(0), fOffset(0),
  fPileUpFile(0), fIndex(0), fBuffer(0),
  fOutputXDR(0), fIndexXDR(0), fBufferXDR(0)
{
  stringstream message;

  fIndex = new char[kIndexSize*8];
  fBuffer = new char[kBufferSize*kRecordSize*4];
  fOutputXDR = new XDR;
  fIndexXDR = new XDR;
  fBufferXDR = new XDR;
  xdrmem_create(fIndexXDR, fIndex, kIndexSize*8, XDR_ENCODE);
  xdrmem_create(fBufferXDR, fBuffer, kBufferSize*kRecordSize*4, XDR_ENCODE);

  fPileUpFile = fopen(fileName, "w+");

  if(fPileUpFile == NULL)
  {
    message << "can't open pile-up file " << fileName;
    throw runtime_error(message.str());
  }

  xdrstdio_create(fOutputXDR, fPileUpFile, XDR_ENCODE);
}

//------------------------------------------------------------------------------

DelphesPileUpWriter::~DelphesPileUpWriter()
{
  xdr_destroy(fOutputXDR);
  if(fPileUpFile) fclose(fPileUpFile);
  xdr_destroy(fBufferXDR);
  xdr_destroy(fIndexXDR);
  if(fBufferXDR) delete fBufferXDR;
  if(fIndexXDR) delete fIndexXDR;
  if(fOutputXDR) delete fOutputXDR;
  if(fBuffer) delete[] fBuffer;
  if(fIndex) delete[] fIndex;
}

//------------------------------------------------------------------------------

void DelphesPileUpWriter::WriteParticle(int pid,
  float x, float y, float z, float t,
  float px, float py, float pz, float e)
{
  if(fEntrySize >= kBufferSize)
  {
    throw runtime_error("too many particles in pile-up event");
  }

  xdr_int(fBufferXDR, &pid);
  xdr_float(fBufferXDR, &x);
  xdr_float(fBufferXDR, &y);
  xdr_float(fBufferXDR, &z);
  xdr_float(fBufferXDR, &t);
  xdr_float(fBufferXDR, &px);
  xdr_float(fBufferXDR, &py);
  xdr_float(fBufferXDR, &pz);
  xdr_float(fBufferXDR, &e);

  ++fEntrySize;
}

//------------------------------------------------------------------------------

void DelphesPileUpWriter::WriteEntry()
{
  if(fEntries >= kIndexSize)
  {
    throw runtime_error("too many pile-up events");
  }

  xdr_int(fOutputXDR, &fEntrySize);
  xdr_opaque(fOutputXDR, fBuffer, fEntrySize*kRecordSize*4);

  xdr_hyper(fIndexXDR, &fOffset);
  fOffset += fEntrySize*kRecordSize*4 + 4;

  xdr_setpos(fBufferXDR, 0);
  fEntrySize = 0;
        
  ++fEntries;
}

//------------------------------------------------------------------------------

void DelphesPileUpWriter::WriteIndex()
{
  xdr_opaque(fOutputXDR, fIndex, fEntries*8);
  xdr_hyper(fOutputXDR, &fEntries);
}

//------------------------------------------------------------------------------
