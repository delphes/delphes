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

//------------------------------------------------------------------------------

#ifndef RecoPuFilter_h
#define RecoPuFilter_h

/** \class RecoPuFilter
 *
 *  Removes particles with RecoPU flag = true.
 *  Input collection needs to pass by TrackPileUpSubtractor first)
 *
 *  \author M. Selvaggi
 *
 */

#include "classes/DelphesModule.h"
#include <vector>

class TIterator;
class TObjArray;

class RecoPuFilter: public DelphesModule
{
public:
  RecoPuFilter();
  ~RecoPuFilter();

  void Init();
  void Process();
  void Finish();

private:
  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(RecoPuFilter, 1)
};

#endif
