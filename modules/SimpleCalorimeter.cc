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

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

SimpleCalorimeter::SimpleCalorimeter() :
  fResolutionFormula(0),
  fItParticleInputArray(0), fItTrackInputArray(0)
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
  set<Double_t>::iterator itPhiBin;
  vector<Double_t> *phiBins;

  // read eta and phi bins
  param = GetParam("EtaPhiBins");
  size = param.GetSize();
  fBinMap.clear();
  fEtaBins.clear();
  fPhiBins.clear();
  for(i = 0; i < size / 2; ++i)
  {
    paramEtaBins = param[i * 2];
    sizeEtaBins = paramEtaBins.GetSize();
    paramPhiBins = param[i * 2 + 1];
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
    phiBins = new vector<double>(itEtaBin->second.size());
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

  for(i = 0; i < size / 2; ++i)
  {
    paramFractions = param[i * 2 + 1];
    fraction = paramFractions[0].GetDouble();
    fFractionMap[param[i * 2].GetInt()] = fraction;
  }

  // read min E value for towers to be saved
  fEnergyMin = GetDouble("EnergyMin", 0.0);

  fEnergySignificanceMin = GetDouble("EnergySignificanceMin", 0.0);

  // flag that says if current calo is Ecal of Hcal (will then fill correct values of Eem and Ehad)
  fIsEcal = GetBool("IsEcal", false);

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

  fEFlowTrackOutputArray = ExportArray(GetString("EFlowTrackOutputArray", "eflowTracks"));
  fEFlowTowerOutputArray = ExportArray(GetString("EFlowTowerOutputArray", "eflowTowers"));
}

//------------------------------------------------------------------------------

void SimpleCalorimeter::Finish()
{
  vector<vector<Double_t> *>::iterator itPhiBin;
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
  Double_t sigma;
  Double_t energyGuess;

  Int_t pdgCode;

  TFractionMap::iterator itFractionMap;

  vector<Double_t>::iterator itEtaBin;
  vector<Double_t>::iterator itPhiBin;
  vector<Double_t> *phiBins;

  vector<Long64_t>::iterator itTowerHits;

  DelphesFactory *factory = GetFactory();
  fTowerHits.clear();
  fTowerFractions.clear();
  fTrackFractions.clear();

  // loop over all particles
  fItParticleInputArray->Reset();
  number = -1;
  fTowerRmax=0.;
  while((particle = static_cast<Candidate *>(fItParticleInputArray->Next())))
  {
    const TLorentzVector &particlePosition = particle->Position;
    ++number;

    // compute maximum radius (needed in FinalizeTower to assess whether barrel or endcap tower)
    if (particlePosition.Perp() > fTowerRmax)
      fTowerRmax=particlePosition.Perp();

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
  while((track = static_cast<Candidate *>(fItTrackInputArray->Next())))
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
    number = (towerHit)&0x0000000000FFFFFFLL;
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
      fTowerEta = 0.5 * (fEtaBins[etaBin - 1] + fEtaBins[etaBin]);
      fTowerPhi = 0.5 * ((*phiBins)[phiBin - 1] + (*phiBins)[phiBin]);

      fTowerEdges[0] = fEtaBins[etaBin - 1];
      fTowerEdges[1] = fEtaBins[etaBin];
      fTowerEdges[2] = (*phiBins)[phiBin - 1];
      fTowerEdges[3] = (*phiBins)[phiBin];

      fTowerEnergy = 0.0;

      fTrackEnergy = 0.0;
      fTrackSigma = 0.0;

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

      track = static_cast<Candidate *>(fTrackInputArray->At(number));
      momentum = track->Momentum;
      position = track->Position;

      energy = momentum.E() * fTrackFractions[number];

      fTrackTime += TMath::Sqrt(energy) * position.T();
      fTrackTimeWeight += TMath::Sqrt(energy);

      if(fTrackFractions[number] > 1.0E-9)
      {

        // compute total charged energy
        fTrackEnergy += energy;
        sigma = fResolutionFormula->Eval(0.0, fTowerEta, 0.0, momentum.E());
        if(sigma / momentum.E() < track->TrackResolution)
          energyGuess = energy;
        else
          energyGuess = momentum.E();

        fTrackSigma += ((track->TrackResolution) * energyGuess) * ((track->TrackResolution) * energyGuess);
        fTowerTrackArray->Add(track);
      }

      else
      {
        fEFlowTrackOutputArray->Add(track);
      }

      continue;
    }

    // check for photon and electron hits in current tower
    if(flags & 2) ++fTowerPhotonHits;

    particle = static_cast<Candidate *>(fParticleInputArray->At(number));
    momentum = particle->Momentum;
    position = particle->Position;

    // fill current tower
    energy = momentum.E() * fTowerFractions[number];

    fTowerEnergy += energy;

    fTowerTime += energy * energy * position.T(); //sigma_t ~ 1/E
    fTowerTimeWeight += energy * energy;

    fTower->AddCandidate(particle);
    fTower->Position = position;

  }

  // finalize last tower
  FinalizeTower();
}

//------------------------------------------------------------------------------

void SimpleCalorimeter::FinalizeTower()
{
  Candidate *tower, *track, *mother;
  Double_t energy, neutralEnergy, pt, eta, phi, r;
  Double_t sigma, neutralSigma;
  Double_t time;

  Double_t weightTrack, weightCalo, bestEnergyEstimate, rescaleFactor;

  TLorentzVector momentum;
  TFractionMap::iterator itFractionMap;

  if(!fTower) return;

  sigma = fResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerEnergy);

  energy = LogNormal(fTowerEnergy, sigma);

  time = (fTowerTimeWeight < 1.0E-09) ? 0.0 : fTowerTime / fTowerTimeWeight;

  sigma = fResolutionFormula->Eval(0.0, fTowerEta, 0.0, energy);

  if(energy < fEnergyMin || energy < fEnergySignificanceMin * sigma) energy = 0.0;

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

  // endcap
  if (TMath::Abs(fTower->Position.Pt() - fTowerRmax) > 1.e-06 && TMath::Abs(eta) > 0.){
    r = fTower->Position.Z()/TMath::SinH(eta);
  }
  // barrel
  else {
    r = fTower->Position.Pt();
  }

  fTower->Position.SetPtEtaPhiE(r, eta, phi, time);
  fTower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
  fTower->L = fTower->Position.Vect().Mag();

  fTower->Eem = (!fIsEcal) ? 0 : energy;
  fTower->Ehad = (fIsEcal) ? 0 : energy;
  fTower->Etrk = fTrackEnergy;

  fTower->Edges[0] = fTowerEdges[0];
  fTower->Edges[1] = fTowerEdges[1];
  fTower->Edges[2] = fTowerEdges[2];
  fTower->Edges[3] = fTowerEdges[3];

  // fill SimpleCalorimeter towers
  if(energy > 0.0) fTowerOutputArray->Add(fTower);

  // e-flow candidates

  //compute neutral excess

  fTrackSigma = TMath::Sqrt(fTrackSigma);
  neutralEnergy = max((energy - fTrackEnergy), 0.0);

  //compute sigma_trk total
  neutralSigma = neutralEnergy / TMath::Sqrt(fTrackSigma * fTrackSigma + sigma * sigma);

  // if neutral excess is significant, simply create neutral Eflow tower and clone each track into eflowtrack
  if(neutralEnergy > fEnergyMin && neutralSigma > fEnergySignificanceMin)
  {
    // create new photon tower
    tower = static_cast<Candidate *>(fTower->Clone());
    pt = neutralEnergy / TMath::CosH(eta);

    tower->Eem = (!fIsEcal) ? 0 : neutralEnergy;
    tower->Ehad = (fIsEcal) ? 0 : neutralEnergy;
    tower->PID = (fIsEcal) ? 22 : 0;

    tower->Momentum.SetPtEtaPhiE(pt, eta, phi, neutralEnergy);
    fEFlowTowerOutputArray->Add(tower);

    fItTowerTrackArray->Reset();
    while((track = static_cast<Candidate *>(fItTowerTrackArray->Next())))
    {
      mother = track;
      track = static_cast<Candidate *>(track->Clone());
      track->AddCandidate(mother);

      fEFlowTrackOutputArray->Add(track);
    }
  }

  // if neutral excess is not significant, rescale eflow tracks, such that the total charged equals the best measurement given by the calorimeter and tracking
  else if(fTrackEnergy > 0.0)
  {
    weightTrack = (fTrackSigma > 0.0) ? 1 / (fTrackSigma * fTrackSigma) : 0.0;
    weightCalo = (sigma > 0.0) ? 1 / (sigma * sigma) : 0.0;

    bestEnergyEstimate = (weightTrack * fTrackEnergy + weightCalo * energy) / (weightTrack + weightCalo);
    rescaleFactor = bestEnergyEstimate / fTrackEnergy;

    fItTowerTrackArray->Reset();
    while((track = static_cast<Candidate *>(fItTowerTrackArray->Next())))
    {
      mother = track;
      track = static_cast<Candidate *>(track->Clone());
      track->AddCandidate(mother);
      track->Momentum.SetPtEtaPhiM(track->Momentum.Pt()*rescaleFactor, track->Momentum.Eta(), track->Momentum.Phi(), track->Momentum.M());
      fEFlowTrackOutputArray->Add(track);
    }
  }
}

//------------------------------------------------------------------------------

Double_t SimpleCalorimeter::LogNormal(Double_t mean, Double_t sigma)
{
  Double_t a, b;

  if(mean > 0.0)
  {
    b = TMath::Sqrt(TMath::Log((1.0 + (sigma * sigma) / (mean * mean))));
    a = TMath::Log(mean) - 0.5 * b * b;

    return TMath::Exp(a + b * gRandom->Gaus(0.0, 1.0));
  }
  else
  {
    return 0.0;
  }
}
