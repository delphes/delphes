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

#ifndef DecayFilter_h
#define DecayFilter_h

/** \class DecayFilter
 *
 *  This module randomly generates decays along the particle trajectory length 
 *  according to actual particle decay length, taking into account for the boost
 *  and using ROOT TDatabasePDG as a source for the particle lifetime.
 *
 *  This module is to to be used after a PropagateParticle step or a similar module
 *  that calculates and store a trajectory length.
 *
 *  Particles that decay are not added to the OutputArray.
 *
 *  \author R. Preghenella - INFN, Bologna
 *
 */

#include "classes/DelphesModule.h"

#include <deque>

class TObjArray;

class DecayFilter: public DelphesModule
{
public:
  DecayFilter();
  ~DecayFilter();

  void Init();
  void Process();
  void Finish();

private:

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(DecayFilter, 1)
};

#endif
