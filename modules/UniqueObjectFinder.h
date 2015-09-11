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

#ifndef UniqueObjectFinder_h
#define UniqueObjectFinder_h

/** \class UniqueObjectFinder
 *
 *  Finds uniquely identified photons, electrons, taus and jets.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <vector>
#include <utility>

class TIterator;
class TObjArray;
class Candidate;

class UniqueObjectFinder: public DelphesModule
{
public:

  UniqueObjectFinder();
  ~UniqueObjectFinder();

  void Init();
  void Process();
  void Finish();

private:

  Bool_t Unique(Candidate *candidate, std::vector< std::pair< TIterator *, TObjArray * > >::iterator itInputMap);

  std::vector< std::pair< TIterator *, TObjArray * > > fInputMap; //!

  ClassDef(UniqueObjectFinder, 1)
};

#endif
