
/** \class DelphesPileUpReader
 *
 *  Reads pile-up binary file
 *
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesPileUpReader.h"

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

DelphesPileUpReader::DelphesPileUpReader(const char *fileName) :
  fEntries(0), fEntrySize(0), fCounter(0),
  fPileUpFile(0), fIndex(0), fBuffer(0),
  fInputXDR(0), fIndexXDR(0), fBufferXDR(0)
{
  stringstream message;

  fIndex = new char[kIndexSize*8];
  fBuffer = new char[kBufferSize*kRecordSize*4];
  fInputXDR = new XDR;
  fIndexXDR = new XDR;
  fBufferXDR = new XDR;
  xdrmem_create(fIndexXDR, fIndex, kIndexSize*8, XDR_DECODE);
  xdrmem_create(fBufferXDR, fBuffer, kBufferSize*kRecordSize*4, XDR_DECODE);

  fPileUpFile = fopen(fileName, "r");

  if(fPileUpFile == NULL)
  {
    message << "can't open pile-up file " << fileName;
    throw runtime_error(message.str());
  }

  xdrstdio_create(fInputXDR, fPileUpFile, XDR_DECODE);

  // read number of events
  fseek(fPileUpFile, -8, SEEK_END);
  xdr_hyper(fInputXDR, &fEntries);

  if(fEntries >= kIndexSize)
  {
    message << "too many events in pile-up file " << fileName;
    throw runtime_error(message.str());
  }

  // read index of events
  fseek(fPileUpFile, -8 - 8*fEntries, SEEK_END);
  xdr_opaque(fInputXDR, fIndex, fEntries*8);
}

//------------------------------------------------------------------------------

DelphesPileUpReader::~DelphesPileUpReader()
{
  xdr_destroy(fInputXDR);
  if(fPileUpFile) fclose(fPileUpFile);
  xdr_destroy(fBufferXDR);
  xdr_destroy(fIndexXDR);
  if(fBufferXDR) delete fBufferXDR;
  if(fIndexXDR) delete fIndexXDR;
  if(fInputXDR) delete fInputXDR;
  if(fBuffer) delete[] fBuffer;
  if(fIndex) delete[] fIndex;
}

//------------------------------------------------------------------------------

bool DelphesPileUpReader::ReadParticle(int &pid,
  float &x, float &y, float &z, float &t,
  float &px, float &py, float &pz, float &e)
{
  if(fCounter >= fEntrySize) return false;

  xdr_int(fBufferXDR, &pid);
  xdr_float(fBufferXDR, &x);
  xdr_float(fBufferXDR, &y);
  xdr_float(fBufferXDR, &z);
  xdr_float(fBufferXDR, &t);
  xdr_float(fBufferXDR, &px);
  xdr_float(fBufferXDR, &py);
  xdr_float(fBufferXDR, &pz);
  xdr_float(fBufferXDR, &e);

  ++fCounter;

  return true;
}

//------------------------------------------------------------------------------

bool DelphesPileUpReader::ReadEntry(quad_t entry)
{
  quad_t offset;

  if(entry >= fEntries) return false;

  // read event position
  xdr_setpos(fIndexXDR, 8*entry);
  xdr_hyper(fIndexXDR, &offset);

  // read event
  fseek(fPileUpFile, offset, SEEK_SET);
  xdr_int(fInputXDR, &fEntrySize);

  if(fEntrySize >= kBufferSize)
  {
    throw runtime_error("too many particles in pile-up event");
  }

  xdr_opaque(fInputXDR, fBuffer, fEntrySize*kRecordSize*4);
  xdr_setpos(fBufferXDR, 0);
  fCounter = 0;

  return true;
}

//------------------------------------------------------------------------------
