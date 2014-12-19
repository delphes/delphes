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

#ifndef ParticlePropagator_h
#define ParticlePropagator_h

/** \class ParticlePropagator
 *
 *  Propagates charged and neutral particles
 *  from a given vertex to a cylinder defined by its radius, 
 *  its half-length, centered at (0,0,0) and with its axis
 *  oriented along the z-axis.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TClonesArray;
class TIterator;

class ParticlePropagator: public DelphesModule
{
public:

  ParticlePropagator();
  ~ParticlePropagator();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fRadius, fRadius2, fHalfLength;
  Double_t fBz;

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!
  TObjArray *fChargedHadronOutputArray; //!
  TObjArray *fElectronOutputArray; //!
  TObjArray *fMuonOutputArray; //!

  ClassDef(ParticlePropagator, 1)
};

#endif
