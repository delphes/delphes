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
#include <memory>
#include <set>
#include <vector>

class TObjArray;
class TVector3;
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
  typedef std::map<Long64_t, std::pair<Double_t, Double_t> > TFractionMap; //!
  typedef std::map<Double_t, std::set<Double_t> > TBinMap; //!

  Candidate *fTower = nullptr;
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

  Bool_t fPhotonAngularSmearing;
  Bool_t fPhotonPointing;

  TFractionMap fFractionMap; //!
  TBinMap fBinMap; //!

  std::vector<Double_t> fEtaBins;
  std::vector<std::unique_ptr<std::vector<Double_t> > > fPhiBins;

  std::vector<Long64_t> fTowerHits;

  std::vector<Double_t> fECalTowerFractions;
  std::vector<Double_t> fHCalTowerFractions;

  std::vector<Double_t> fECalTrackFractions;
  std::vector<Double_t> fHCalTrackFractions;

  std::unique_ptr<DelphesFormula> fECalResolutionFormula; //!
  std::unique_ptr<DelphesFormula> fHCalResolutionFormula; //!

  std::unique_ptr<DelphesFormula> fPhotonThetaResolutionFormula; //!
  std::unique_ptr<DelphesFormula> fPhotonPhiResolutionFormula; //!
  std::unique_ptr<DelphesFormula> fPointingThetaResolutionFormula; //!
  std::unique_ptr<DelphesFormula> fPointingPhiResolutionFormula; //!

  std::unique_ptr<TIterator> fItParticleInputArray; //!
  std::unique_ptr<TIterator> fItTrackInputArray; //!

  const TObjArray *fParticleInputArray = nullptr; //!
  const TObjArray *fTrackInputArray = nullptr; //!

  TObjArray *fTowerOutputArray = nullptr; //!
  TObjArray *fPhotonOutputArray = nullptr; //!

  TObjArray *fEFlowTrackOutputArray = nullptr; //!
  TObjArray *fEFlowPhotonOutputArray = nullptr; //!
  TObjArray *fEFlowNeutralHadronOutputArray = nullptr; //!

  std::unique_ptr<TObjArray> fECalTowerTrackArray; //!
  std::unique_ptr<TIterator> fItECalTowerTrackArray; //!

  std::unique_ptr<TObjArray> fHCalTowerTrackArray; //!
  std::unique_ptr<TIterator> fItHCalTowerTrackArray; //!

  std::unique_ptr<TObjArray> fTowerTrackArray; //!
  std::unique_ptr<TIterator> fItTowerTrackArray; //!

  void FinalizeTower();
  Bool_t SmearPhotonDirection(Candidate *candidate, const TVector3 &impact, const TVector3 &flight, Double_t energy, Double_t sigmaE, Double_t &etaOut, Double_t &phiOut);
  Double_t LogNormal(Double_t mean, Double_t sigma);
  Double_t TruncatedGaussian(Double_t mean, Double_t sigma);

  ClassDef(DualReadoutCalorimeter, 1)
};

#endif
