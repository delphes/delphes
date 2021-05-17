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

#ifndef TimeOfFlight_h
#define TimeOfFlight_h

/** \class TimeOfFlight
 *
 *  Calculates Time-Of-Flight
 *
 *  \author Michele Selvaggi - CERN
 *
*/

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;

class TimeOfFlight: public DelphesModule
{
public:
  TimeOfFlight();
  ~TimeOfFlight();

  void Init();
  void Process();
  void Finish();
  void ComputeVertexMomenta();

private:

  Int_t fVertexTimeMode;

  TIterator *fItTrackInputArray; //!
  TIterator *fItVertexInputArray; //!

  const TObjArray *fTrackInputArray; //!
  const TObjArray *fVertexInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(TimeOfFlight, 1)
};

#endif
