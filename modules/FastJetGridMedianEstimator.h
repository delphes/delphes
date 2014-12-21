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

#ifndef FastJetGridMedianEstimator_h
#define FastJetGridMedianEstimator_h


/** \class FastJetGridMedianEstimator
 *
 *  Computes median energy density per event using a fixed grid.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"
#include <vector>

class TObjArray;
class TIterator;

namespace fastjet {
  class GridMedianBackgroundEstimator;
}

class FastJetGridMedianEstimator: public DelphesModule
{
public:

  FastJetGridMedianEstimator();
  ~FastJetGridMedianEstimator();

  void Init();
  void Process();
  void Finish();

private:

  std::vector< fastjet::GridMedianBackgroundEstimator * > fEstimators; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fRhoOutputArray; //!

  ClassDef(FastJetGridMedianEstimator, 1)
};

#endif
