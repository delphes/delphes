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

/** \class DelphesXDRWriter
 *
 *  Writes XDR
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesXDRWriter.h"

#include <array>

//------------------------------------------------------------------------------

void DelphesXDRWriter::SetBuffer(void *buffer)
{
  fBuffer = reinterpret_cast<uint8_t *>(buffer);
  fOffset = 0;
}

//------------------------------------------------------------------------------

void DelphesXDRWriter::WriteRaw(void *value, int size)
{
  int rndup = size % 4;
  if(rndup > 0) rndup = 4 - rndup;
  if(fFile)
    fwrite(value, 1, size + rndup, fFile);
}

//------------------------------------------------------------------------------

void DelphesXDRWriter::WriteValue(void *value, int size)
{
  const uint8_t *src = reinterpret_cast<uint8_t *>(value);
  if(fBuffer)
  {
    for(int i = 0; i < size; ++i) fBuffer[fOffset + i] = src[size - 1 - i];
    fOffset += size;
  }
  else if(fFile)
  {
    std::array<uint8_t, 8> buffer;
    for(int i = 0; i < size; ++i) buffer[i] = src[size - 1 - i];
    WriteRaw(buffer.data(), size);
  }
  else
    return;
}

//------------------------------------------------------------------------------
