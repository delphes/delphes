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

#include <array>

#include <string.h>

//------------------------------------------------------------------------------

void DelphesXDRReader::SetBuffer(void *buffer)
{
  fBuffer = reinterpret_cast<uint8_t *>(buffer);
  fOffset = 0;
}

//------------------------------------------------------------------------------

void DelphesXDRReader::ReadRaw(void *value, int size)
{
  int rndup = size % 4;
  if(rndup > 0)
    rndup = 4 - rndup;
  if(fFile)
    fread(value, 1, size + rndup, fFile);
}

//------------------------------------------------------------------------------

void DelphesXDRReader::ReadValue(void *value, int size)
{
  uint8_t *dst = reinterpret_cast<uint8_t *>(value);
  if(fBuffer)
  {
    fOffset += size;
    for(int i = 0; i < size; ++i) dst[i] = fBuffer[fOffset - 1 - i];
  }
  else if(fFile)
  {
    std::array<uint8_t, 8> buffer;
    ReadRaw(buffer.data(), size);
    for(int i = 0; i < size; ++i) dst[i] = buffer.at(size - 1 - i);
  }
}

//------------------------------------------------------------------------------

void DelphesXDRReader::ReadString(void *value, int maxSize)
{
  int32_t size;
  if(ReadValue(&size, 4); size < 0 || size > maxSize) size = maxSize;

  if(fBuffer)
  {
    memcpy(value, fBuffer + fOffset, size);
    fOffset += size;
  }
  else if(fFile)
    ReadRaw(value, size);
}

//------------------------------------------------------------------------------
