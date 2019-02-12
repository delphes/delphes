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

#ifndef DelphesXDRWriter_h
#define DelphesXDRWriter_h

/** \class DelphesXDRWriter
 *
 *  Writes XDR
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include <stdint.h>
#include <stdio.h>

class DelphesXDRWriter
{
public:
  DelphesXDRWriter();

  void SetFile(FILE *file);
  void SetBuffer(void *buffer);
  void SetOffset(int offset);

  void WriteRaw(void *value, int size);
  void WriteValue(void *value, int size);

private:
  FILE *fFile;
  uint8_t *fBuffer;

  int fOffset;
};

#endif // DelphesXDRWriter_h
