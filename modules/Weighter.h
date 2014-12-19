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

#ifndef Weighter_h
#define Weighter_h

/** \class Weighter
 *
 *  Apply a weight depending on PDG code.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <set>
#include <map>

class TObjArray;

class Weighter: public DelphesModule
{
public:

  Weighter();
  ~Weighter();

  void Init();
  void Process();
  void Finish();

private:

#if !defined(__CINT__) && !defined(__CLING__)
  struct TIndexStruct
  {
    Int_t codes[4];
    bool operator< (const TIndexStruct &value) const;
  };

  std::set<Int_t> fWeightSet, fCodeSet;
  std::map<TIndexStruct, Double_t> fWeightMap;
#endif

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(Weighter, 1)
};

#endif
