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

#ifndef Calorimeter_h
#define Calorimeter_h

/** \class Calorimeter
 *
 *  Fills calorimeter towers, performs calorimeter resolution smearing,
 *  and creates energy flow objects (tracks, photons, and neutral hadrons).
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <map>
#include <set>
#include <vector>

class TObjArray;
class DelphesFormula;
class Candidate;

class Calorimeter: public DelphesModule
{
public:

  Calorimeter();
  ~Calorimeter();

  void Init();
  void Process();
  void Finish();

private:

  typedef std::map< Long64_t, std::pair< Double_t, Double_t > > TFractionMap; //!
  typedef std::map< Double_t, std::set< Double_t > > TBinMap; //!

  Candidate *fTower;
  Double_t fTowerEta, fTowerPhi, fTowerEdges[4];
  Double_t fECalTowerEnergy, fHCalTowerEnergy;
  Double_t fECalTrackEnergy[2], fHCalTrackEnergy[2];

  Double_t fTimingEnergyMin;
  Bool_t fElectronsFromTrack;

  Int_t fTowerTrackHits, fTowerPhotonHits;

  Double_t fECalEnergyMin;
  Double_t fHCalEnergyMin;

  Double_t fECalEnergySignificanceMin;
  Double_t fHCalEnergySignificanceMin;

  Bool_t fSmearTowerCenter;

  TFractionMap fFractionMap; //!
  TBinMap fBinMap; //!

  std::vector < Double_t > fEtaBins;
  std::vector < std::vector < Double_t >* > fPhiBins;

  std::vector < Long64_t > fTowerHits;

  std::vector < Double_t > fECalTowerFractions;
  std::vector < Double_t > fHCalTowerFractions;

  std::vector < Double_t > fECalTrackFractions;
  std::vector < Double_t > fHCalTrackFractions;

  DelphesFormula *fECalResolutionFormula; //!
  DelphesFormula *fHCalResolutionFormula; //!

  TIterator *fItParticleInputArray; //!
  TIterator *fItTrackInputArray; //!

  const TObjArray *fParticleInputArray; //!
  const TObjArray *fTrackInputArray; //!

  TObjArray *fTowerOutputArray; //!
  TObjArray *fPhotonOutputArray; //!

  TObjArray *fEFlowTrackOutputArray; //!
  TObjArray *fEFlowPhotonOutputArray; //!
  TObjArray *fEFlowNeutralHadronOutputArray; //!

  TObjArray *fECalTowerTrackArray[2]; //!
  TIterator *fItECalTowerTrackArray[2]; //!

  TObjArray *fHCalTowerTrackArray[2]; //!
  TIterator *fItHCalTowerTrackArray[2]; //!

  void FinalizeTower();
  Double_t LogNormal(Double_t mean, Double_t sigma);

  ClassDef(Calorimeter, 1)
};

#endif
