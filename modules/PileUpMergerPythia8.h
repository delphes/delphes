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

#ifndef PileUpMergerPythia8_h
#define PileUpMergerPythia8_h

/** \class PileUpMergerPythia8
 *
 *  Merges particles from pile-up sample into event
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TObjArray;
class DelphesTF2;

namespace Pythia8
{
 class Pythia;
};

class PileUpMergerPythia8: public DelphesModule
{
public:

  PileUpMergerPythia8();
  ~PileUpMergerPythia8();

  void Init();
  void Process();
  void Finish();

private:

  Int_t fPileUpDistribution;
  Double_t fMeanPileUp;

  Double_t fZVertexSpread;
  Double_t fTVertexSpread;

  Double_t fInputBeamSpotX;
  Double_t fInputBeamSpotY;
  Double_t fOutputBeamSpotX;
  Double_t fOutputBeamSpotY;

  Double_t fPTMin;

  DelphesTF2 *fFunction; //!

  Pythia8::Pythia *fPythia; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fParticleOutputArray; //!
  TObjArray *fVertexOutputArray; //!

  ClassDef(PileUpMergerPythia8, 1)
};

#endif
