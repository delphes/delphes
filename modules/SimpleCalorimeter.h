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

#ifndef SimpleCalorimeter_h
#define SimpleCalorimeter_h

/** \class SimpleCalorimeter
 *
 *  Fills SimpleCalorimeter towers, performs SimpleCalorimeter resolution smearing,
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

class SimpleCalorimeter: public DelphesModule
{
public:

  SimpleCalorimeter();
  ~SimpleCalorimeter();

  void Init();
  void Process();
  void Finish();

private:

  typedef std::map< Long64_t, Double_t > TFractionMap; //!
  typedef std::map< Double_t, std::set< Double_t > > TBinMap; //!

  Candidate *fTower;
  Double_t fTowerEta, fTowerPhi, fTowerEdges[4];
  Double_t fTowerEnergy;
  Double_t fTrackEnergy;

  Double_t fTowerTime;
  Double_t fTrackTime;

  Double_t fTowerTimeWeight;
  Double_t fTrackTimeWeight;

  Int_t fTowerTrackHits, fTowerPhotonHits;

  Double_t fEnergyMin;

  Double_t fEnergySignificanceMin;

  Bool_t fSmearTowerCenter;

  TFractionMap fFractionMap; //!
  TBinMap fBinMap; //!

  std::vector < Double_t > fEtaBins;
  std::vector < std::vector < Double_t >* > fPhiBins;

  std::vector < Long64_t > fTowerHits;

  std::vector < Double_t > fTowerFractions;

  std::vector < Double_t > fTrackFractions;

  DelphesFormula *fResolutionFormula; //!

  TIterator *fItParticleInputArray; //!
  TIterator *fItTrackInputArray; //!

  const TObjArray *fParticleInputArray; //!
  const TObjArray *fTrackInputArray; //!

  TObjArray *fTowerOutputArray; //!

  TObjArray *fEFlowTowerOutputArray; //!

  TObjArray *fTowerTrackArray; //!
  TIterator *fItTowerTrackArray; //!

  void FinalizeTower();
  Double_t LogNormal(Double_t mean, Double_t sigma);

  ClassDef(SimpleCalorimeter, 1)
};

#endif
