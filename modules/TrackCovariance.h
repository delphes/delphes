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

#ifndef TrackCovariance_h
#define TrackCovariance_h

/** \class TrackCovariance
 *
 *  Smears track parameters according to appropriate covariance matrix.
 *
 *  \authors P. Demin - UCLouvain, Louvain-la-Neuve
 *           M. Selvaggi - CERN
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;

class SolGeom;
class SolGridCov;
class AcceptanceClx;
class DelphesFormula;

class TrackCovariance : public DelphesModule
{
public:
  TrackCovariance();
  ~TrackCovariance();

  void Init();
  void Process();
  void Finish();

private:
  Double_t fBz;
  Int_t fNMinHits;

  const std::unique_ptr<DelphesFormula> fElectronScaleFactor;
  const std::unique_ptr<DelphesFormula> fMuonScaleFactor;
  const std::unique_ptr<DelphesFormula> fChargedHadronScaleFactor;

  const std::unique_ptr<SolGeom> fGeometry;
  const std::unique_ptr<SolGridCov> fCovariance;

  AcceptanceClx *fAcx{nullptr};

  TIterator *fItInputArray{nullptr}; //!

  const TObjArray *fInputArray{nullptr}; //!

  TObjArray *fOutputArray{nullptr}; //!

  ClassDef(TrackCovariance, 1)
};

#endif
