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

#ifndef PileUpMerger_h
#define PileUpMerger_h

/** \class PileUpMerger
 *
 *  Merges particles from pile-up sample into event
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TObjArray;
class DelphesPileUpReader;
class DelphesTF2;

class PileUpMerger: public DelphesModule
{
public:

  PileUpMerger();
  ~PileUpMerger();

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

  DelphesTF2 *fFunction; //!

  DelphesPileUpReader *fReader; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fParticleOutputArray; //!
  TObjArray *fVertexOutputArray; //!

  ClassDef(PileUpMerger, 1)
};

#endif
