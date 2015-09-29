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

#ifndef PhotonConversions_h
#define PhotonConversions_h

/** \class PhotonConversions
 *
 *  Converts photons into e+ e- pairs according to material ditribution in the detector.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TClonesArray;
class TIterator;
class DelphesCylindricalFormula;
class TF1;

class PhotonConversions: public DelphesModule
{
public:

  PhotonConversions();
  ~PhotonConversions();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fRadius, fRadius2, fHalfLength;
  Double_t fEtaMin, fEtaMax;

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  DelphesCylindricalFormula *fConversionMap; //!

  TF1 *fDecayXsec; //!

  Double_t fStep;

  ClassDef(PhotonConversions, 1)
};

#endif
