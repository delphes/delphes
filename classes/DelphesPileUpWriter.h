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

#ifndef DelphesPileUpWriter_h
#define DelphesPileUpWriter_h

/** \class DelphesPileUpWriter
 *
 *  Writes pile-up binary file
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include <stdio.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

class DelphesPileUpWriter
{
public:

  DelphesPileUpWriter(const char *fileName);

  ~DelphesPileUpWriter();

  void WriteParticle(int pid,
    float x, float y, float z, float t,
    float px, float py, float pz, float e);

  void WriteEntry();

  void WriteIndex();

private:

  quad_t fEntries;
  int fEntrySize;
  quad_t fOffset;

  FILE *fPileUpFile;
  char *fIndex;
  char *fBuffer;

  XDR *fOutputXDR;
  XDR *fIndexXDR;
  XDR *fBufferXDR;
};

#endif // DelphesPileUpWriter_h


