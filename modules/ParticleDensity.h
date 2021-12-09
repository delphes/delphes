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

#ifndef ParticleDensity_h
#define ParticleDensity_h

/** \class ParticleDensity
 *
 *  This module calculates the particle multiplicity density in eta-phi bins.
 *  It then assigns the value to the candidates according to the candidate eta.
 *
 *  \author R. Preghenella - INFN, Bologna
 *
 */

#include "classes/DelphesModule.h"

#include <deque>

class TObjArray;
class TH2F;

class ParticleDensity: public DelphesModule
{
public:
  ParticleDensity();
  ~ParticleDensity();

  void Init();
  void Process();
  void Finish();

private:

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  Bool_t fUseMomentumVector; // !
  TH2F *fHisto; //!
  
  ClassDef(ParticleDensity, 1)
};

#endif
