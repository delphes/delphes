/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2020  Universite catholique de Louvain (UCLouvain), Belgium
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

#ifndef ClusterCounting_h
#define ClusterCounting_h

/** \class ClusterCounting
 *
 *  Counts ionisation clusters of energy loss in drift chambers
 *
 *  \authors F. Bedeschi - INFN
 *           M. Selvaggi - CERN
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;

class TrkUtil;

class ClusterCounting: public DelphesModule
{
public:
  ClusterCounting();
  ~ClusterCounting();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fRmin;
	Double_t fRmax;
	Double_t fZmin;
	Double_t fZmax;
  Double_t fBz;

  Int_t fGasOption;

  TrkUtil *fTrackUtil;

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(ClusterCounting, 1)
};

#endif
