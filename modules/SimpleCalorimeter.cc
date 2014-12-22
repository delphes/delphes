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


/** \class SimpleCalorimeter
 *
 *  Fills SimpleCalorimeter towers, performs SimpleCalorimeter resolution smearing,
 *  and creates energy flow objects (tracks, photons, and neutral hadrons).
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/SimpleCalorimeter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

SimpleCalorimeter::SimpleCalorimeter() :
  fResolutionFormula(0),
  fItParticleInputArray(0), fItTrackInputArray(0),
  fTowerTrackArray(0), fItTowerTrackArray(0)
{
  fResolutionFormula = new DelphesFormula;

  fTowerTrackArray = new TObjArray;
  fItTowerTrackArray = fTowerTrackArray->MakeIterator();
}

//------------------------------------------------------------------------------

SimpleCalorimeter::~SimpleCalorimeter()
{
  if(fResolutionFormula) delete fResolutionFormula;

  if(fTowerTrackArray) delete fTowerTrackArray;
  if(fItTowerTrackArray) delete fItTowerTrackArray;
}

//------------------------------------------------------------------------------

void SimpleCalorimeter::Init()
{
  ExRootConfParam param, paramEtaBins, paramPhiBins, paramFractions;
  Long_t i, j, k, size, sizeEtaBins, sizePhiBins;
  Double_t fraction;
  TBinMap::iterator itEtaBin;
  set< Double_t >::iterator itPhiBin;
  vector< Double_t > *phiBins;

  // read eta and phi bins
  param = GetParam("EtaPhiBins");
  size = param.GetSize();
  fBinMap.clear();
  fEtaBins.clear();
  fPhiBins.clear();
  for(i = 0; i < size/2; ++i)
  {
    paramEtaBins = param[i*2];
    sizeEtaBins = paramEtaBins.GetSize();
    paramPhiBins = param[i*2 + 1];
    sizePhiBins = paramPhiBins.GetSize();

    for(j = 0; j < sizeEtaBins; ++j)
    {
      for(k = 0; k < sizePhiBins; ++k)
      {
        fBinMap[paramEtaBins[j].GetDouble()].insert(paramPhiBins[k].GetDouble());
      }
    }
  }

  // for better performance we transform map of sets to parallel vectors:
  // vector< double > and vector< vector< double >* >
  for(itEtaBin = fBinMap.begin(); itEtaBin != fBinMap.end(); ++itEtaBin)
  {
    fEtaBins.push_back(itEtaBin->first);
    phiBins = new vector< double >(itEtaBin->second.size());
    fPhiBins.push_back(phiBins);
    phiBins->clear();
    for(itPhiBin = itEtaBin->second.begin(); itPhiBin != itEtaBin->second.end(); ++itPhiBin)
    {
      phiBins->push_back(*itPhiBin);
    }
  }

  // read energy fractions for different particles
  param = GetParam("EnergyFraction");
  size = param.GetSize();

  // set default energy fractions values
  fFractionMap.clear();
  fFractionMap[0] = 1.0;

  for(i = 0; i < size/2; ++i)
  {
    paramFractions = param[i*2 + 1];
    fraction = paramFractions[0].GetDouble();
    fFractionMap[param[i*2].GetInt()] = fraction;
  }

/*
  TFractionMap::iterator itFractionMap;
  for(itFractionMap = fFractionMap.begin(); itFractionMap != fFractionMap.end(); ++itFractionMap)
  {
    cout << itFractionMap->first << "   " << itFractionMap->second.first  << "   " << itFractionMap->second.second << endl;
  }
*/

  // read min E value for towers to be saved
  fEnergyMin = GetDouble("EnergyMin", 0.0);

  fEnergySignificanceMin = GetDouble("EnergySignificanceMin", 0.0);

  // switch on or off the dithering of the center of calorimeter towers
  fSmearTowerCenter = GetBool("SmearTowerCenter", true);

  // read resolution formulas
  fResolutionFormula->Compile(GetString("ResolutionFormula", "0"));

  // import array with output from other modules
  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "ParticlePropagator/particles"));
  fItParticleInputArray = fParticleInputArray->MakeIterator();

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "ParticlePropagator/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  // create output arrays
  fTowerOutputArray = ExportArray(GetString("TowerOutputArray", "towers"));
  fEFlowTowerOutputArray = ExportArray(GetString("EFlowTowerOutputArray", "eflowTowers"));
}

//------------------------------------------------------------------------------

void SimpleCalorimeter::Finish()
{
  vector< vector< Double_t >* >::iterator itPhiBin;
  if(fItParticleInputArray) delete fItParticleInputArray;
  if(fItTrackInputArray) delete fItTrackInputArray;
  for(itPhiBin = fPhiBins.begin(); itPhiBin != fPhiBins.end(); ++itPhiBin)
  {
    delete *itPhiBin;
  }
}

//------------------------------------------------------------------------------

void SimpleCalorimeter::Process()
{
  Candidate *particle, *track;
  TLorentzVector position, momentum;
  Short_t etaBin, phiBin, flags;
  Int_t number;
  Long64_t towerHit, towerEtaPhi, hitEtaPhi;
  Double_t fraction;
  Double_t energy;
  Int_t pdgCode;

  TFractionMap::iterator itFractionMap;

  vector< Double_t >::iterator itEtaBin;
  vector< Double_t >::iterator itPhiBin;
  vector< Double_t > *phiBins;

  vector< Long64_t >::iterator itTowerHits;

  DelphesFactory *factory = GetFactory();
  fTowerHits.clear();
  fTowerFractions.clear();
  fTrackFractions.clear();

  // loop over all particles
  fItParticleInputArray->Reset();
  number = -1;
  while((particle = static_cast<Candidate*>(fItParticleInputArray->Next())))
  {
    const TLorentzVector &particlePosition = particle->Position;
    ++number;

    pdgCode = TMath::Abs(particle->PID);

    itFractionMap = fFractionMap.find(pdgCode);
    if(itFractionMap == fFractionMap.end())
    {
      itFractionMap = fFractionMap.find(0);
    }

    fraction = itFractionMap->second;
    fTowerFractions.push_back(fraction);

    if(fraction < 1.0E-9) continue;

    // find eta bin [1, fEtaBins.size - 1]
    itEtaBin = lower_bound(fEtaBins.begin(), fEtaBins.end(), particlePosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    etaBin = distance(fEtaBins.begin(), itEtaBin);

    // phi bins for given eta bin
    phiBins = fPhiBins[etaBin];

    // find phi bin [1, phiBins.size - 1]
    itPhiBin = lower_bound(phiBins->begin(), phiBins->end(), particlePosition.Phi());
    if(itPhiBin == phiBins->begin() || itPhiBin == phiBins->end()) continue;
    phiBin = distance(phiBins->begin(), itPhiBin);

    flags = 0;
    flags |= (pdgCode == 11 || pdgCode == 22) << 1;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for particle number}
    towerHit = (Long64_t(etaBin) << 48) | (Long64_t(phiBin) << 32) | (Long64_t(flags) << 24) | Long64_t(number);

    fTowerHits.push_back(towerHit);
  }

  // loop over all tracks
  fItTrackInputArray->Reset();
  number = -1;
  while((track = static_cast<Candidate*>(fItTrackInputArray->Next())))
  {
    const TLorentzVector &trackPosition = track->Position;
    ++number;

    pdgCode = TMath::Abs(track->PID);

    itFractionMap = fFractionMap.find(pdgCode);
    if(itFractionMap == fFractionMap.end())
    {
      itFractionMap = fFractionMap.find(0);
    }

    fraction = itFractionMap->second;

    fTrackFractions.push_back(fraction);

    // find eta bin [1, fEtaBins.size - 1]
    itEtaBin = lower_bound(fEtaBins.begin(), fEtaBins.end(), trackPosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    etaBin = distance(fEtaBins.begin(), itEtaBin);

    // phi bins for given eta bin
    phiBins = fPhiBins[etaBin];

    // find phi bin [1, phiBins.size - 1]
    itPhiBin = lower_bound(phiBins->begin(), phiBins->end(), trackPosition.Phi());
    if(itPhiBin == phiBins->begin() || itPhiBin == phiBins->end()) continue;
    phiBin = distance(phiBins->begin(), itPhiBin);

    flags = 1;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for track number}
    towerHit = (Long64_t(etaBin) << 48) | (Long64_t(phiBin) << 32) | (Long64_t(flags) << 24) | Long64_t(number);

    fTowerHits.push_back(towerHit);
  }

  // all hits are sorted first by eta bin number, then by phi bin number,
  // then by flags and then by particle or track number
  sort(fTowerHits.begin(), fTowerHits.end());

  // loop over all hits
  towerEtaPhi = 0;
  fTower = 0;
  for(itTowerHits = fTowerHits.begin(); itTowerHits != fTowerHits.end(); ++itTowerHits)
  {
    towerHit = (*itTowerHits);
    flags = (towerHit >> 24) & 0x00000000000000FFLL;
    number = (towerHit) & 0x0000000000FFFFFFLL;
    hitEtaPhi = towerHit >> 32;

    if(towerEtaPhi != hitEtaPhi)
    {
      // switch to next tower
      towerEtaPhi = hitEtaPhi;

      // finalize previous tower
      FinalizeTower();

      // create new tower
      fTower = factory->NewCandidate();

      phiBin = (towerHit >> 32) & 0x000000000000FFFFLL;
      etaBin = (towerHit >> 48) & 0x000000000000FFFFLL;

      // phi bins for given eta bin
      phiBins = fPhiBins[etaBin];

      // calculate eta and phi of the tower's center
      fTowerEta = 0.5*(fEtaBins[etaBin - 1] + fEtaBins[etaBin]);
      fTowerPhi = 0.5*((*phiBins)[phiBin - 1] + (*phiBins)[phiBin]);

      fTowerEdges[0] = fEtaBins[etaBin - 1];
      fTowerEdges[1] = fEtaBins[etaBin];
      fTowerEdges[2] = (*phiBins)[phiBin - 1];
      fTowerEdges[3] = (*phiBins)[phiBin];

      fTowerEnergy = 0.0;
      fTrackEnergy = 0.0;

      fTowerTime = 0.0;
      fTrackTime = 0.0;

      fTowerTimeWeight = 0.0;

      fTowerTrackHits = 0;
      fTowerPhotonHits = 0;

      fTowerTrackArray->Clear();
    }

    // check for track hits
    if(flags & 1)
    {
      ++fTowerTrackHits;

      track = static_cast<Candidate*>(fTrackInputArray->At(number));
      momentum = track->Momentum;
      position = track->Position;

      energy = momentum.E() * fTrackFractions[number];

      fTrackEnergy += energy;

      fTrackTime += TMath::Sqrt(energy)*position.T();
      fTrackTimeWeight += TMath::Sqrt(energy);

      fTowerTrackArray->Add(track);

      continue;
    }

    // check for photon and electron hits in current tower
    if(flags & 2) ++fTowerPhotonHits;

    particle = static_cast<Candidate*>(fParticleInputArray->At(number));
    momentum = particle->Momentum;
    position = particle->Position;

    // fill current tower
    energy = momentum.E() * fTowerFractions[number];

    fTowerEnergy += energy;

    fTowerTime += TMath::Sqrt(energy)*position.T();
    fTowerTimeWeight += TMath::Sqrt(energy);

    fTower->AddCandidate(particle);
  }

  // finalize last tower
  FinalizeTower();
}

//------------------------------------------------------------------------------

void SimpleCalorimeter::FinalizeTower()
{
  Candidate *tower;
  Double_t energy, pt, eta, phi;
  Double_t sigma;
  Double_t time;

  if(!fTower) return;

  sigma = fResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerEnergy);

  energy = LogNormal(fTowerEnergy, sigma);

  time = (fTowerTimeWeight < 1.0E-09 ) ? 0.0 : fTowerTime/fTowerTimeWeight;

  sigma = fResolutionFormula->Eval(0.0, fTowerEta, 0.0, energy);

  if(energy < fEnergyMin || energy < fEnergySignificanceMin*sigma) energy = 0.0;

  if(fSmearTowerCenter)
  {
    eta = gRandom->Uniform(fTowerEdges[0], fTowerEdges[1]);
    phi = gRandom->Uniform(fTowerEdges[2], fTowerEdges[3]);
  }
  else
  {
    eta = fTowerEta;
    phi = fTowerPhi;
  }

  pt = energy / TMath::CosH(eta);

  fTower->Position.SetPtEtaPhiE(1.0, eta, phi, time);
  fTower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);

  fTower->Edges[0] = fTowerEdges[0];
  fTower->Edges[1] = fTowerEdges[1];
  fTower->Edges[2] = fTowerEdges[2];
  fTower->Edges[3] = fTowerEdges[3];

  // fill SimpleCalorimeter towers
  if(energy > 0.0) fTowerOutputArray->Add(fTower);

  // fill energy flow candidates
  energy -= fTrackEnergy;

  sigma = fResolutionFormula->Eval(0.0, fTowerEta, 0.0, energy);

  if(energy < fEnergyMin || energy < fEnergySignificanceMin*sigma) energy = 0.0;

  // save energy excess as an energy flow tower
  if(energy > 0.0)
  {
    // create new photon tower
    tower = static_cast<Candidate*>(fTower->Clone());
    pt = energy / TMath::CosH(eta);

    tower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
    fEFlowTowerOutputArray->Add(tower);
  }
}

//------------------------------------------------------------------------------

Double_t SimpleCalorimeter::LogNormal(Double_t mean, Double_t sigma)
{
  Double_t a, b;

  if(mean > 0.0)
  {
    b = TMath::Sqrt(TMath::Log((1.0 + (sigma*sigma)/(mean*mean))));
    a = TMath::Log(mean) - 0.5*b*b;

    return TMath::Exp(a + b*gRandom->Gaus(0.0, 1.0));
  }
  else
  {
    return 0.0;
  }
}
