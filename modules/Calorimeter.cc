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


/** \class Calorimeter
 *
 *  Fills calorimeter towers, performs calorimeter resolution smearing,
 *  and creates energy flow objects (tracks, photons, and neutral hadrons).
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Calorimeter.h"

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

Calorimeter::Calorimeter() :
  fECalResolutionFormula(0), fHCalResolutionFormula(0),
  fItParticleInputArray(0), fItTrackInputArray(0)
{
  Int_t i;

  fECalResolutionFormula = new DelphesFormula;
  fHCalResolutionFormula = new DelphesFormula;

  for(i = 0; i < 2; ++i)
  {
    fECalTowerTrackArray[i] = new TObjArray;
    fItECalTowerTrackArray[i] = fECalTowerTrackArray[i]->MakeIterator();

    fHCalTowerTrackArray[i] = new TObjArray;
    fItHCalTowerTrackArray[i] = fHCalTowerTrackArray[i]->MakeIterator();
  }
}

//------------------------------------------------------------------------------

Calorimeter::~Calorimeter()
{
  Int_t i;

  if(fECalResolutionFormula) delete fECalResolutionFormula;
  if(fHCalResolutionFormula) delete fHCalResolutionFormula;

  for(i = 0; i < 2; ++i)
  {
    if(fECalTowerTrackArray[i]) delete fECalTowerTrackArray[i];
    if(fItECalTowerTrackArray[i]) delete fItECalTowerTrackArray[i];

    if(fHCalTowerTrackArray[i]) delete fHCalTowerTrackArray[i];
    if(fItHCalTowerTrackArray[i]) delete fItHCalTowerTrackArray[i];
  }
}

//------------------------------------------------------------------------------

void Calorimeter::Init()
{
  ExRootConfParam param, paramEtaBins, paramPhiBins, paramFractions;
  Long_t i, j, k, size, sizeEtaBins, sizePhiBins;
  Double_t ecalFraction, hcalFraction;
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
  fFractionMap[0] = make_pair(0.0, 1.0);

  for(i = 0; i < size/2; ++i)
  {
    paramFractions = param[i*2 + 1];

    ecalFraction = paramFractions[0].GetDouble();
    hcalFraction = paramFractions[1].GetDouble();

    fFractionMap[param[i*2].GetInt()] = make_pair(ecalFraction, hcalFraction);
  }

  // read min E value for timing measurement in ECAL
  fTimingEnergyMin = GetDouble("TimingEnergyMin",4.);
  // For timing
  // So far this flag needs to be false
  // Curved extrapolation not supported
  fElectronsFromTrack = false;

  // read min E value for towers to be saved
  fECalEnergyMin = GetDouble("ECalEnergyMin", 0.0);
  fHCalEnergyMin = GetDouble("HCalEnergyMin", 0.0);

  fECalEnergySignificanceMin = GetDouble("ECalEnergySignificanceMin", 0.0);
  fHCalEnergySignificanceMin = GetDouble("HCalEnergySignificanceMin", 0.0);

  // switch on or off the dithering of the center of calorimeter towers
  fSmearTowerCenter = GetBool("SmearTowerCenter", true);

  // read resolution formulas
  fECalResolutionFormula->Compile(GetString("ECalResolutionFormula", "0"));
  fHCalResolutionFormula->Compile(GetString("HCalResolutionFormula", "0"));

  // import array with output from other modules
  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "ParticlePropagator/particles"));
  fItParticleInputArray = fParticleInputArray->MakeIterator();

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "ParticlePropagator/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  // create output arrays
  fTowerOutputArray = ExportArray(GetString("TowerOutputArray", "towers"));
  fPhotonOutputArray = ExportArray(GetString("PhotonOutputArray", "photons"));

  fEFlowTrackOutputArray = ExportArray(GetString("EFlowTrackOutputArray", "eflowTracks"));
  fEFlowPhotonOutputArray = ExportArray(GetString("EFlowPhotonOutputArray", "eflowPhotons"));
  fEFlowNeutralHadronOutputArray = ExportArray(GetString("EFlowNeutralHadronOutputArray", "eflowNeutralHadrons"));
}

//------------------------------------------------------------------------------

void Calorimeter::Finish()
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

void Calorimeter::Process()
{
  Candidate *particle, *track;
  TLorentzVector position, momentum;
  Short_t etaBin, phiBin, flags;
  Int_t number;
  Long64_t towerHit, towerEtaPhi, hitEtaPhi;
  Double_t ecalFraction, hcalFraction;
  Double_t ecalEnergy, hcalEnergy;
  Double_t ecalSigma, hcalSigma;
  Int_t pdgCode;

  TFractionMap::iterator itFractionMap;

  vector< Double_t >::iterator itEtaBin;
  vector< Double_t >::iterator itPhiBin;
  vector< Double_t > *phiBins;

  vector< Long64_t >::iterator itTowerHits;

  DelphesFactory *factory = GetFactory();
  fTowerHits.clear();
  fECalTowerFractions.clear();
  fHCalTowerFractions.clear();
  fECalTrackFractions.clear();
  fHCalTrackFractions.clear();

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

    ecalFraction = itFractionMap->second.first;
    hcalFraction = itFractionMap->second.second;

    fECalTowerFractions.push_back(ecalFraction);
    fHCalTowerFractions.push_back(hcalFraction);

    if(ecalFraction < 1.0E-9 && hcalFraction < 1.0E-9) continue;

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

    ecalFraction = itFractionMap->second.first;
    hcalFraction = itFractionMap->second.second;

    fECalTrackFractions.push_back(ecalFraction);
    fHCalTrackFractions.push_back(hcalFraction);

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

      fECalTowerEnergy = 0.0;
      fHCalTowerEnergy = 0.0;

      fECalTrackEnergy[0] = 0.0;
      fECalTrackEnergy[1] = 0.0;

      fHCalTrackEnergy[0] = 0.0;
      fHCalTrackEnergy[1] = 0.0;

      fTowerTrackHits = 0;
      fTowerPhotonHits = 0;

      fECalTowerTrackArray[0]->Clear();
      fECalTowerTrackArray[1]->Clear();

      fHCalTowerTrackArray[0]->Clear();
      fHCalTowerTrackArray[1]->Clear();
    }

    // check for track hits
    if(flags & 1)
    {
      ++fTowerTrackHits;

      track = static_cast<Candidate*>(fTrackInputArray->At(number));
      momentum = track->Momentum;
      position = track->Position;

      ecalEnergy = momentum.E() * fECalTrackFractions[number];
      hcalEnergy = momentum.E() * fHCalTrackFractions[number];

      if(ecalEnergy > fTimingEnergyMin && fTower)
      {
        if(fElectronsFromTrack)
        {
          fTower->ECalEnergyTimePairs.push_back(make_pair<Float_t, Float_t>(ecalEnergy, track->Position.T()));
        }
      }

      if(fECalTrackFractions[number] > 1.0E-9 && fHCalTrackFractions[number] < 1.0E-9)
      {
        ecalSigma = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, momentum.E());
        if(ecalSigma/momentum.E() < track->TrackResolution)
        {
          fECalTrackEnergy[0] += ecalEnergy;
          fECalTowerTrackArray[0]->Add(track);
        }
        else
        {
          fECalTrackEnergy[1] += ecalEnergy;
          fECalTowerTrackArray[1]->Add(track);
        }
      }
      else if(fECalTrackFractions[number] < 1.0E-9 && fHCalTrackFractions[number] > 1.0E-9)
      {
        hcalSigma = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, momentum.E());
        if(hcalSigma/momentum.E() < track->TrackResolution)
        {
          fHCalTrackEnergy[0] += hcalEnergy;
          fHCalTowerTrackArray[0]->Add(track);
        }
        else
        {
          fHCalTrackEnergy[1] += hcalEnergy;
          fHCalTowerTrackArray[1]->Add(track);
        }
      }
      else if(fECalTrackFractions[number] < 1.0E-9 && fHCalTrackFractions[number] < 1.0E-9)
      {
        fEFlowTrackOutputArray->Add(track);
      }

      continue;
    }

    // check for photon and electron hits in current tower
    if(flags & 2) ++fTowerPhotonHits;

    particle = static_cast<Candidate*>(fParticleInputArray->At(number));
    momentum = particle->Momentum;
    position = particle->Position;

    // fill current tower
    ecalEnergy = momentum.E() * fECalTowerFractions[number];
    hcalEnergy = momentum.E() * fHCalTowerFractions[number];

    fECalTowerEnergy += ecalEnergy;
    fHCalTowerEnergy += hcalEnergy;

    if(ecalEnergy > fTimingEnergyMin && fTower)
    {
      if (abs(particle->PID) != 11 || !fElectronsFromTrack)
      {
        fTower->ECalEnergyTimePairs.push_back(make_pair<Float_t, Float_t>(ecalEnergy, particle->Position.T()));
      }
    }

    fTower->AddCandidate(particle);
  }

  // finalize last tower
  FinalizeTower();
}

//------------------------------------------------------------------------------

void Calorimeter::FinalizeTower()
{
  Candidate *track, *tower, *mother;
  Double_t energy, pt, eta, phi;
  Double_t ecalEnergy, hcalEnergy;
  Double_t ecalSigma, hcalSigma;

  TLorentzVector momentum;
  TFractionMap::iterator itFractionMap;

  Float_t weight, sumWeightedTime, sumWeight;

  if(!fTower) return;


  ecalSigma = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fECalTowerEnergy);
  hcalSigma = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fHCalTowerEnergy);

  ecalEnergy = LogNormal(fECalTowerEnergy, ecalSigma);
  hcalEnergy = LogNormal(fHCalTowerEnergy, hcalSigma);

  ecalSigma = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, ecalEnergy);
  hcalSigma = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, hcalEnergy);

  if(ecalEnergy < fECalEnergyMin || ecalEnergy < fECalEnergySignificanceMin*ecalSigma) ecalEnergy = 0.0;
  if(hcalEnergy < fHCalEnergyMin || hcalEnergy < fHCalEnergySignificanceMin*hcalSigma) hcalEnergy = 0.0;

  energy = ecalEnergy + hcalEnergy;

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

  // Time calculation for tower
  fTower->NTimeHits = 0;
  sumWeightedTime = 0.0;
  sumWeight = 0.0;

  for(size_t i = 0; i < fTower->ECalEnergyTimePairs.size(); ++i)
  {
    weight = TMath::Sqrt(fTower->ECalEnergyTimePairs[i].first);
    sumWeightedTime += weight * fTower->ECalEnergyTimePairs[i].second;
    sumWeight += weight;
    fTower->NTimeHits++;
  }

  if(sumWeight > 0.0)
  {
    fTower->Position.SetPtEtaPhiE(1.0, eta, phi, sumWeightedTime/sumWeight);
  }
  else
  {
    fTower->Position.SetPtEtaPhiE(1.0, eta, phi, 999999.9);
  }


  fTower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
  fTower->Eem = ecalEnergy;
  fTower->Ehad = hcalEnergy;

  fTower->Edges[0] = fTowerEdges[0];
  fTower->Edges[1] = fTowerEdges[1];
  fTower->Edges[2] = fTowerEdges[2];
  fTower->Edges[3] = fTowerEdges[3];

  if(energy > 0.0)
  {
    if(fTowerPhotonHits > 0 && fTowerTrackHits == 0)
    {
      fPhotonOutputArray->Add(fTower);
    }

    fTowerOutputArray->Add(fTower);
  }

  // fill energy flow candidates

  ecalEnergy -= fECalTrackEnergy[1];
  hcalEnergy -= fHCalTrackEnergy[1];

  fItECalTowerTrackArray[0]->Reset();
  while((track = static_cast<Candidate*>(fItECalTowerTrackArray[0]->Next())))
  {
    mother = track;
    track = static_cast<Candidate*>(track->Clone());
    track->AddCandidate(mother);

    track->Momentum *= ecalEnergy/fECalTrackEnergy[0];

    fEFlowTrackOutputArray->Add(track);
  }

  fItECalTowerTrackArray[1]->Reset();
  while((track = static_cast<Candidate*>(fItECalTowerTrackArray[1]->Next())))
  {
    mother = track;
    track = static_cast<Candidate*>(track->Clone());
    track->AddCandidate(mother);

    fEFlowTrackOutputArray->Add(track);
  }

  fItHCalTowerTrackArray[0]->Reset();
  while((track = static_cast<Candidate*>(fItHCalTowerTrackArray[0]->Next())))
  {
    mother = track;
    track = static_cast<Candidate*>(track->Clone());
    track->AddCandidate(mother);

    track->Momentum *= hcalEnergy/fHCalTrackEnergy[0];

    fEFlowTrackOutputArray->Add(track);
  }

  fItHCalTowerTrackArray[1]->Reset();
  while((track = static_cast<Candidate*>(fItHCalTowerTrackArray[1]->Next())))
  {
    mother = track;
    track = static_cast<Candidate*>(track->Clone());
    track->AddCandidate(mother);

    fEFlowTrackOutputArray->Add(track);
  }

  if(fECalTowerTrackArray[0]->GetEntriesFast() > 0) ecalEnergy = 0.0;
  if(fHCalTowerTrackArray[0]->GetEntriesFast() > 0) hcalEnergy = 0.0;

  ecalSigma = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, ecalEnergy);
  hcalSigma = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, hcalEnergy);

  if(ecalEnergy < fECalEnergyMin || ecalEnergy < fECalEnergySignificanceMin*ecalSigma) ecalEnergy = 0.0;
  if(hcalEnergy < fHCalEnergyMin || hcalEnergy < fHCalEnergySignificanceMin*hcalSigma) hcalEnergy = 0.0;

  energy = ecalEnergy + hcalEnergy;

  if(ecalEnergy > 0.0)
  {
    // create new photon tower
    tower = static_cast<Candidate*>(fTower->Clone());

    pt = ecalEnergy / TMath::CosH(eta);

    tower->Momentum.SetPtEtaPhiE(pt, eta, phi, ecalEnergy);
    tower->Eem = ecalEnergy;
    tower->Ehad = 0.0;

    fEFlowPhotonOutputArray->Add(tower);
  }
  if(hcalEnergy > 0.0)
  {
    // create new neutral hadron tower
    tower = static_cast<Candidate*>(fTower->Clone());

    pt = hcalEnergy / TMath::CosH(eta);

    tower->Momentum.SetPtEtaPhiE(pt, eta, phi, hcalEnergy);
    tower->Eem = 0.0;
    tower->Ehad = hcalEnergy;

    fEFlowNeutralHadronOutputArray->Add(tower);
  }
}

//------------------------------------------------------------------------------

Double_t Calorimeter::LogNormal(Double_t mean, Double_t sigma)
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
