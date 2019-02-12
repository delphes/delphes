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

#ifndef PhotonID_h
#define PhotonID_h

/** \class PhotonID
 *
 *  Applies complex photon Id. Reconstructed photon candidtes are first separated into matched and non-matched to gen particles. 
 *  Non-matched pass the "fake" efficiency. Matched photons get further splitted into isolated and non-isolated (user can choose criterion for isolation)
 *  Isolated photons pass the "prompt" efficiency while the non-isolated pass the "non-prompt" efficiency
 *
 *  \author M. Selvaggi - CERN
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;
class DelphesFormula;
class Candidate;

class PhotonID: public DelphesModule
{
public:
  PhotonID();
  ~PhotonID();

  void Init();
  void Process();
  void Finish();

private:
  DelphesFormula *fPromptFormula;
  DelphesFormula *fNonPromptFormula;
  DelphesFormula *fFakeFormula;

  // import input arrays
  const TObjArray *fInputPhotonArray;
  TIterator *fItInputPhotonArray;

  // use filtered collection for speed
  const TObjArray *fInputGenArray;
  TIterator *fItInputGenArray;

  Double_t fPTMin;
  Double_t fRelIsoMax;

  TObjArray *fOutputArray; //!

  Bool_t isFake(const Candidate *obj);

  ClassDef(PhotonID, 1)
};

#endif
