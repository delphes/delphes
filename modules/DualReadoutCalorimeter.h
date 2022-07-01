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

#ifndef DualReadoutCalorimeter_h
#define DualReadoutCalorimeter_h

/** \class DualReadoutCalorimeter
 *
 *  Fills DualReadoutCalorimeter towers, performs DualReadoutCalorimeter resolution smearing,
 *  and creates energy flow objects (tracks, photons, and neutral hadrons).
 *
 *  \author M. Selvaggi - CERN
 *
 */

#include "classes/DelphesModule.h"

#include <map>
#include <set>
#include <vector>

class TObjArray;
class DelphesFormula;
class Candidate;

class DualReadoutCalorimeter: public DelphesModule
{
public:

  DualReadoutCalorimeter();
  ~DualReadoutCalorimeter();

  void Init();
  void Process();
  void Finish();

private:

  typedef std::map< Long64_t, std::pair< Double_t, Double_t > > TFractionMap; //!
  typedef std::map< Double_t, std::set< Double_t > > TBinMap; //!

  Candidate *fTower;
  Double_t fTowerEta, fTowerPhi, fTowerEdges[4];
  Double_t fECalTowerEnergy, fHCalTowerEnergy;
  Double_t fECalTrackEnergy, fHCalTrackEnergy;
  Double_t fTrackEnergy;
  Double_t fTowerRmax;

  Double_t fTimingEnergyMin;
  Bool_t fElectronsFromTrack;

  Int_t fTowerTrackHits, fTowerPhotonHits;

  Double_t fECalEnergyMin;
  Double_t fHCalEnergyMin;
  Double_t fEnergyMin;

  Double_t fECalMinSignificance;
  Double_t fHCalMinSignificance;

  Double_t fEnergySignificanceMin;

  Double_t fECalTrackSigma;
  Double_t fHCalTrackSigma;
  Double_t fTrackSigma;

  Double_t fTowerTime;
  Double_t fTowerTimeWeight;

  Bool_t fSmearTowerCenter;
  Bool_t fSmearLogNormal;

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

  TObjArray *fECalTowerTrackArray; //!
  TIterator *fItECalTowerTrackArray; //!

  TObjArray *fHCalTowerTrackArray; //!
  TIterator *fItHCalTowerTrackArray; //!

  TObjArray *fTowerTrackArray; //!
  TIterator *fItTowerTrackArray; //!

  void FinalizeTower();
  Double_t LogNormal(Double_t mean, Double_t sigma);
  Double_t TruncatedGaussian(Double_t mean, Double_t sigma);

  ClassDef(DualReadoutCalorimeter, 1)
};

#endif
