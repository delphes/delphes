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

#ifndef TrackPileUpSubtractor_h
#define TrackPileUpSubtractor_h

/** \class TrackPileUpSubtractor
 *
 *  Subtract pile-up contribution from tracks.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <map>

class TIterator;
class TObjArray;

class TrackPileUpSubtractor: public DelphesModule
{
public:

  TrackPileUpSubtractor();
  ~TrackPileUpSubtractor();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fZVertexResolution;

  Double_t fPTMin; 

  std::map< TIterator *, TObjArray * > fInputMap; //!

  ClassDef(TrackPileUpSubtractor, 1)

  TIterator *fItVertexInputArray; //!

  const TObjArray *fVertexInputArray; //!

};

#endif
