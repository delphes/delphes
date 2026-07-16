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

#ifndef BeamSpotSmearing_h
#define BeamSpotSmearing_h

/** \class BeamSpotSmearing
 *  Applies per-event Gaussian smearing to the beam spot and the bunch time spread.
 *
 *  \author Michele Selvaggi - CERN
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;

class BeamSpotSmearing: public DelphesModule
{
public:
  BeamSpotSmearing();
  ~BeamSpotSmearing();

  void Init();
  void Process();
  void Finish();

private:
  Double_t fSigmaX; //! beam spot size in x, stored in mm
  Double_t fSigmaY; //! beam spot size in y, stored in mm
  Double_t fSigmaZ; //! beam spot size in z, stored in mm
  Double_t fSigmaT; //! bunch time spread, stored in mm (c*t)

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  ClassDef(BeamSpotSmearing, 1)
};

#endif
