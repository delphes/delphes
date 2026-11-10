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

#ifndef EventCounterFilter_h
#define EventCounterFilter_h

/** \class EventCounterFilter
 *
 *  Selects events by object multiplicity.
 *
 *  \author Sihyun Jeon (Boston University)
 *
 */

#include "classes/DelphesModule.h"

#include <string>
#include <vector>

class TObjArray;
class DelphesTreeWriter;

class EventCounterFilter: public DelphesModule
{
public:
  EventCounterFilter();
  ~EventCounterFilter();

  void Init();
  void Process();
  void Finish();

private:
  enum Op
  {
    kEq,
    kNe,
    kGt,
    kGe,
    kLt,
    kLe
  };

  std::vector<TObjArray *> fInputArrays; //!
  std::vector<Op> fOps; //!
  std::vector<Int_t> fCounts; //!
  std::vector<std::string> fNames; //!

  Long64_t fTotal; //!
  Long64_t fPassed; //!

  DelphesTreeWriter *fWriter; //!

  ClassDef(EventCounterFilter, 1)
};

#endif
