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

#ifndef DelphesPileUpReader_h
#define DelphesPileUpReader_h

/** \class DelphesPileUpReader
 *
 *  Reads pile-up binary file
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include <stdio.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

class DelphesPileUpReader
{
public:

  DelphesPileUpReader(const char *fileName);

  ~DelphesPileUpReader();

  bool ReadParticle(int &pid,
    float &x, float &y, float &z, float &t,
    float &px, float &py, float &pz, float &e);

  bool ReadEntry(quad_t entry);

  quad_t GetEntries() const { return fEntries; }

private:

  quad_t fEntries;

  int fEntrySize;
  int fCounter;

  FILE *fPileUpFile;
  char *fIndex;
  char *fBuffer;

  XDR *fInputXDR;
  XDR *fIndexXDR;
  XDR *fBufferXDR;
};

#endif // DelphesPileUpReader_h


