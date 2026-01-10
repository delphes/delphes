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
 *  \author A. Chattopadhyay - UPRM (Added insensitive calorimeter bins functionality)
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
  Long_t i, j, k, size, sizeEtaBins, sizePhiBins;
  Double_t fraction;
  TBinMap::iterator itEtaBin;
  set<Double_t>::iterator itPhiBin;
  vector<Double_t> *phiBins;

  // read eta and phi bins
  auto param = GetParam("EtaPhiBins");
  size = param->GetSize();
  fBinMap.clear();
  fEtaBins.clear();
  fPhiBins.clear();
  for(i = 0; i < size / 2; ++i)
  {
    const auto paramEtaBins = (*param)[i * 2];
    sizeEtaBins = paramEtaBins->GetSize();
    const auto paramPhiBins = (*param)[i * 2 + 1];
    sizePhiBins = paramPhiBins->GetSize();

    for(j = 0; j < sizeEtaBins; ++j)
    {
      for(k = 0; k < sizePhiBins; ++k)
      {
        fBinMap[(*paramEtaBins)[j]->GetDouble()].insert((*paramPhiBins)[k]->GetDouble());
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

  // for blind calorimeter: read insensitive bins (convert centers -> indices)
  fInsensitiveBinSet.clear();

  // Get insensitive bin parameters
  const auto blindParam = GetParam("InsensitiveEtaPhiBins");
  Long_t nBlind = blindParam->GetSize();

  // Loop over blind eta-phi pairs
  for(Long_t ib = 0; ib < nBlind; ++ib)
  {
    const auto pairParam = (*blindParam)[ib];
    if(pairParam->GetSize() < 2) continue;

    double etaEdge = (*pairParam)[0]->GetDouble();
    double phiEdge = (*pairParam)[1]->GetDouble();

    // Find closest eta bin index (using bin centers)
    auto itEta = std::min_element(fEtaBins.begin(), fEtaBins.end(),
                                  [etaEdge](double a, double b) {
                                      return std::abs(a - etaEdge) < std::abs(b - etaEdge);
                                  });
    if(itEta == fEtaBins.end()) continue;
    Short_t etaBinIndex = static_cast<Short_t>(std::distance(fEtaBins.begin(), itEta));

    // Find closest phi bin index (using bin centers) 
    std::vector<double>* phiVec = fPhiBins[etaBinIndex];
    if (!phiVec || phiVec->empty()) continue;

    auto itPhi = std::min_element(phiVec->begin(), phiVec->end(),
                                  [phiEdge](double a, double b) {
                                      return std::abs(a - phiEdge) < std::abs(b - phiEdge);
                                  });
    if(itPhi == phiVec->end()) continue;
    Short_t phiBinIndex = static_cast<Short_t>(std::distance(phiVec->begin(), itPhi));

    // Insert bin into insensitive set 
    fInsensitiveBinSet.insert(std::make_pair(etaBinIndex, phiBinIndex));
  }

  // Prompt to let user know how many insensitive bins were saved
  std::cout << "SimpleBlindCalorimeter: insensitive bins saved = " << fInsensitiveBinSet.size() << std::endl;

  // read energy fractions for different particles
  param = GetParam("EnergyFraction");
  size = param->GetSize();

  // set default energy fractions values
  fFractionMap.clear();
  fFractionMap[0] = 1.0;

  for(i = 0; i < size / 2; ++i)
  {
    const auto paramFractions = (*param)[i * 2 + 1];
    fraction = (*paramFractions)[0]->GetDouble();
    fFractionMap[(*param)[i * 2]->GetInt()] = fraction;
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
    

    // skip insensitive calorimeter bins entirely for particles
    if (IsTowerInsensitive(etaBin, phiBin))  {
      fTower = nullptr;
   } 
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
    // skip insensitive calo bins entirely for particles
    if (IsTowerInsensitive(etaBin, phiBin)) {
      fTower = nullptr;
    } 
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
          if(fTower) FinalizeTower();

          // create new tower
          fTower = factory->NewCandidate();
          phiBin = (towerHit >> 32) & 0x000000000000FFFFLL;
          etaBin = (towerHit >> 48) & 0x000000000000FFFFLL;

          //mark fTower nullptr
          if (IsTowerInsensitive(etaBin, phiBin))  {
            fTower = nullptr;  // insensitive tower: no creation
            // do NOT continue here! preserve hit loop for ordering
          }

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
          fNeutralEnergy = 0.0;
          fTrackSigma = 0.0;

          fTowerEnergyFromPU = 0.0;
          fTrackEnergyFromPU = 0.0;
          fNeutralEnergyFromPU = 0.0;

          fTowerTime = 0.0;
          fTrackTime = 0.0;

          fTowerTimeWeight = 0.0;

          fTowerTrackHits = 0;
          fTowerPhotonHits = 0;

          fTowerTrackArray->Clear();
      }

      // We already handled it above; fTower == nullptr if insensitive
      // This prevents skipping hits that are already stored in fTowerHits


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

            // compute energy contribution from PU tracks
            if (track->IsPU) fTrackEnergyFromPU += energy;
          
            sigma = fResolutionFormula->Eval(0.0, fTowerEta, 0.0, momentum.E());
            if(sigma / momentum.E() < track->TrackResolution)
                energyGuess = energy;
            else
                energyGuess = momentum.E();

            fTrackSigma += ((track->TrackResolution) * energyGuess) * ((track->TrackResolution) * energyGuess);
            if(fTower) fTowerTrackArray->Add(track);  // add only if tower exists
          }
          else
          {
            if(fTower) fEFlowTrackOutputArray->Add(track);
          }
        continue;
      }
      else
      {

        particle = static_cast<Candidate *>(fParticleInputArray->At(number));
        momentum = particle->Momentum;

        energy = momentum.E() * fTrackFractions[number];

        if(fTrackFractions[number] > 1.0E-9)
        {
            // compute total neutral energy
            fNeutralEnergy += energy;
            if (particle->IsPU) fNeutralEnergyFromPU += energy;
        }     
      }

    // fill current tower
    energy = momentum.E() * fTowerFractions[number];

    // check for photon and electron hits in current tower
    if(flags & 2) ++fTowerPhotonHits;

    particle = static_cast<Candidate *>(fParticleInputArray->At(number));
    momentum = particle->Momentum;
    position = particle->Position;


    // fill current tower
    energy = momentum.E() * fTowerFractions[number];
    if(fTower)  // add only if tower exists
    {
        fTowerEnergy += energy;
        fTowerTime += energy * energy * position.T(); //sigma_t ~ 1/E
        fTowerTimeWeight += energy * energy;
        fTower->AddCandidate(particle);
        fTower->Position = position;
        if (particle->IsPU) fTowerEnergyFromPU += energy; 
    }

  } 

  // finalize last tower (only if it was sensitive)
  if(fTower) FinalizeTower();
    
}

//------------------------------------------------------------------------------

void SimpleCalorimeter::FinalizeTower()
{
  Candidate *tower, *track, *mother;
  Double_t energy, neutralEnergy, pt, eta, phi, r;
  Double_t sigma, neutralSigma;
  Double_t time;

  Double_t weightTrack, weightCalo, bestEnergyEstimate, rescaleFactor;

  Double_t hardEnergyFraction;

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

  hardEnergyFraction = (fTowerEnergy > 1.0E-09) ? (fTowerEnergy - fTowerEnergyFromPU) / fTowerEnergy : 1.0;
  
  fTower->BetaStar = hardEnergyFraction;

  // fill  SimpleCalorimeter towers
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

    // compute hard energy fraction for tower
    hardEnergyFraction = (fNeutralEnergy > 1.0E-09) ? (fNeutralEnergy - fNeutralEnergyFromPU) / (fNeutralEnergy) : 1.0;

    if (hardEnergyFraction < 0.5) tower->IsPU = 1;
    else tower->IsPU = 0;

    tower->BetaStar = hardEnergyFraction;
    
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
