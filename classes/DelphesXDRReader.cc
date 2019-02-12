/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2018  Universite catholique de Louvain (UCL), Belgium
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

/** \class DelphesXDRReader
 *
 *  Reads XDR
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesXDRReader.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

//------------------------------------------------------------------------------

DelphesXDRReader::DelphesXDRReader() :
  fFile(0), fBuffer(0), fOffset(0)
{
}

//------------------------------------------------------------------------------

void DelphesXDRReader::SetFile(FILE *file)
{
  fFile = file;
}

//------------------------------------------------------------------------------

void DelphesXDRReader::SetBuffer(void *buffer)
{
  fBuffer = (uint8_t *)buffer;
  fOffset = 0;
}

//------------------------------------------------------------------------------

void DelphesXDRReader::SetOffset(int offset)
{
  fOffset = offset;
}

//------------------------------------------------------------------------------

void DelphesXDRReader::ReadRaw(void *value, int size)
{
  int rndup;

  rndup = size % 4;
  if(rndup > 0)
  {
    rndup = 4 - rndup;
  }

  if(fFile)
  {
    fread(value, 1, size + rndup, fFile);
  }
}

//------------------------------------------------------------------------------

void DelphesXDRReader::ReadValue(void *value, int size)
{
  int i;
  uint8_t *dst, buffer[8];

  dst = (uint8_t *)value;

  if(fBuffer)
  {
    fOffset += size;
    for(i = 0; i < size; ++i) dst[i] = fBuffer[fOffset - 1 - i];
  }
  else if(fFile)
  {
    ReadRaw(buffer, size);
    for(i = 0; i < size; ++i) dst[i] = buffer[size - 1 - i];
  }
}

//------------------------------------------------------------------------------

void DelphesXDRReader::ReadString(void *value, int maxSize)
{
  int32_t size;

  ReadValue(&size, 4);

  if(size > maxSize) size = maxSize;

  if(fBuffer)
  {
    memcpy(value, fBuffer + fOffset, size);
    fOffset += size;
  }
  else if(fFile)
  {
    ReadRaw(value, size);
  }
}

//------------------------------------------------------------------------------
