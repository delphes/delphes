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

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

#include <algorithm>
#include <map>
#include <set>
#include <vector>

using namespace std;

class Calorimeter: public DelphesModule
{
public:
  explicit Calorimeter(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fTimingEnergyMin(Steer<double>("TimingEnergyMin", 4.)), // read min E value for timing measurement in ECAL
    fECalEnergyMin(Steer<double>("ECalEnergyMin", 0.0)), // read min E value for towers to be saved
    fHCalEnergyMin(Steer<double>("HCalEnergyMin", 0.0)),
    fECalEnergySignificanceMin(Steer<double>("ECalEnergySignificanceMin", 0.0)),
    fHCalEnergySignificanceMin(Steer<double>("HCalEnergySignificanceMin", 0.0)),
    fSmearTowerCenter(Steer<bool>("SmearTowerCenter", true)), // switch on or off the dithering of the center of calorimeter towers
    fPhiBins(Steer<std::unordered_map<double, std::vector<double> > >("EtaPhiBins")),
    fFractionMap(Steer<TFractionMap>("EnergyFraction")), // ECAL/HCAL energy fractions
    fECalTowerTrackArray(std::make_shared<std::vector<Candidate *> >()),
    fHCalTowerTrackArray(std::make_shared<std::vector<Candidate *> >()),
    fECalResolutionFormula(std::make_unique<DelphesFormula>()),
    fHCalResolutionFormula(std::make_unique<DelphesFormula>())
  {
    for(const std::pair<double, std::vector<double> > etaPhiBins : fPhiBins) // auto would avoid a copy
      fEtaBins.emplace_back(etaPhiBins.first);
    std::sort(fEtaBins.begin(), fEtaBins.end());

    // set default energy fractions values
    fFractionMap[0] = std::make_pair(0., 1.);

    // read resolution formulas
    fECalResolutionFormula->Compile(Steer<std::string>("ECalResolutionFormula", "0"));
    fHCalResolutionFormula->Compile(Steer<std::string>("HCalResolutionFormula", "0"));
  }

  void Init() override
  {
    // import array with output from other modules
    fParticleInputArray = ImportArray(Steer<std::string>("ParticleInputArray", "ParticlePropagator/particles"));
    fTrackInputArray = ImportArray(Steer<std::string>("TrackInputArray", "ParticlePropagator/tracks"));
    // create output arrays
    fTowerOutputArray = ExportArray(Steer<std::string>("TowerOutputArray", "towers"));
    fPhotonOutputArray = ExportArray(Steer<std::string>("PhotonOutputArray", "photons"));
    fEFlowTrackOutputArray = ExportArray(Steer<std::string>("EFlowTrackOutputArray", "eflowTracks"));
    fEFlowPhotonOutputArray = ExportArray(Steer<std::string>("EFlowPhotonOutputArray", "eflowPhotons"));
    fEFlowNeutralHadronOutputArray = ExportArray(Steer<std::string>("EFlowNeutralHadronOutputArray", "eflowNeutralHadrons"));
  }
  void Process() override;

private:
  typedef std::map<Long64_t, std::pair<double, Double_t> > TFractionMap; //!
  typedef std::map<double, std::set<Double_t> > TBinMap; //!

  void FinalizeTower();
  double LogNormal(Double_t mean, Double_t sigma);

  const double fTimingEnergyMin;
  // For timing
  // So far this flag needs to be false
  // Curved extrapolation not supported
  const bool fElectronsFromTrack{false};
  const double fECalEnergyMin;
  const double fHCalEnergyMin;

  const double fECalEnergySignificanceMin;
  const double fHCalEnergySignificanceMin;

  const bool fSmearTowerCenter;

  const std::unordered_map<double, std::vector<double> > fPhiBins;
  std::vector<double> fEtaBins;

  TFractionMap fFractionMap; //!

  // unstored collections
  const CandidatesCollection fECalTowerTrackArray; //!
  const CandidatesCollection fHCalTowerTrackArray; //!

  const std::unique_ptr<DelphesFormula> fECalResolutionFormula; //!
  const std::unique_ptr<DelphesFormula> fHCalResolutionFormula; //!

  double fTowerRmax;
  int fTowerTrackHits, fTowerPhotonHits;

  double fECalTrackSigma;
  double fHCalTrackSigma;

  Candidate *fTower{nullptr};
  double fTowerEta, fTowerPhi, fTowerEdges[4];
  double fECalTowerEnergy, fHCalTowerEnergy;
  double fECalTrackEnergy, fHCalTrackEnergy;

  TBinMap fBinMap; //!

  std::vector<Long64_t> fTowerHits;

  std::vector<double> fECalTowerFractions;
  std::vector<double> fHCalTowerFractions;

  std::vector<double> fECalTrackFractions;
  std::vector<double> fHCalTrackFractions;

  // input collections
  CandidatesCollection fParticleInputArray; //!
  CandidatesCollection fTrackInputArray; //!

  // stored output collections
  CandidatesCollection fTowerOutputArray; //!
  CandidatesCollection fPhotonOutputArray; //!

  CandidatesCollection fEFlowTrackOutputArray; //!
  CandidatesCollection fEFlowPhotonOutputArray; //!
  CandidatesCollection fEFlowNeutralHadronOutputArray; //!
};

//------------------------------------------------------------------------------

void Calorimeter::Process()
{
  fTowerOutputArray->clear();
  fPhotonOutputArray->clear();
  fEFlowTrackOutputArray->clear();
  fEFlowPhotonOutputArray->clear();
  fEFlowNeutralHadronOutputArray->clear();

  Candidate *particle, *track;
  TLorentzVector position, momentum;
  Short_t etaBin, phiBin, flags;
  int number;
  Long64_t towerHit, towerEtaPhi, hitEtaPhi;
  double ecalFraction, hcalFraction;
  double ecalEnergy, hcalEnergy;
  double ecalSigma, hcalSigma;
  double energyGuess;
  int pdgCode;

  TFractionMap::iterator itFractionMap;

  vector<Long64_t>::iterator itTowerHits;

  DelphesFactory *factory = GetFactory();
  fTowerHits.clear();
  fECalTowerFractions.clear();
  fHCalTowerFractions.clear();
  fECalTrackFractions.clear();
  fHCalTrackFractions.clear();

  // loop over all particles
  number = -1;
  fTowerRmax = 0.;
  for(Candidate *const &particle : *fParticleInputArray)
  {
    const TLorentzVector &particlePosition = particle->Position;
    ++number;

    // compute maximum radius (needed in FinalizeTower to assess whether barrel or endcap tower)
    if(particlePosition.Perp() > fTowerRmax)
      fTowerRmax = particlePosition.Perp();

    pdgCode = std::abs(particle->PID);

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
    std::vector<double>::const_iterator itEtaBin = std::lower_bound(fEtaBins.begin(), fEtaBins.end(), particlePosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    etaBin = std::distance(fEtaBins.cbegin(), itEtaBin);

    // phi bins for given eta bin
    const std::vector<double> &phiBins = fPhiBins.at(*itEtaBin);

    // find phi bin [1, phiBins.size - 1]
    std::vector<double>::const_iterator itPhiBin = std::lower_bound(phiBins.begin(), phiBins.end(), particlePosition.Phi());
    if(itPhiBin == phiBins.begin() || itPhiBin == phiBins.end()) continue;
    phiBin = std::distance(phiBins.cbegin(), itPhiBin);

    flags = 0;
    flags |= (pdgCode == 11 || pdgCode == 22) << 1;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for particle number}
    towerHit = (Long64_t(etaBin) << 48) | (Long64_t(phiBin) << 32) | (Long64_t(flags) << 24) | Long64_t(number);

    fTowerHits.push_back(towerHit);
  }

  // loop over all tracks
  number = -1;
  for(Candidate *const &track : *fTrackInputArray)
  {
    const TLorentzVector &trackPosition = track->Position;
    ++number;

    pdgCode = std::abs(track->PID);

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
    std::vector<double>::const_iterator itEtaBin = std::lower_bound(fEtaBins.begin(), fEtaBins.end(), trackPosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    etaBin = distance(fEtaBins.cbegin(), itEtaBin);

    // phi bins for given eta bin
    const std::vector<double> &phiBins = fPhiBins.at(*itEtaBin);

    // find phi bin [1, phiBins.size - 1]
    std::vector<double>::const_iterator itPhiBin = std::lower_bound(phiBins.begin(), phiBins.end(), trackPosition.Phi());
    if(itPhiBin == phiBins.begin() || itPhiBin == phiBins.end()) continue;
    phiBin = std::distance(phiBins.cbegin(), itPhiBin);

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
      const std::vector<double> &phiBins = fPhiBins.at(fEtaBins.at(etaBin));

      // calculate eta and phi of the tower's center
      fTowerEta = 0.5 * (fEtaBins.at(etaBin - 1) + fEtaBins.at(etaBin));
      fTowerPhi = 0.5 * (phiBins.at(phiBin - 1) + phiBins.at(phiBin));

      fTowerEdges[0] = fEtaBins.at(etaBin - 1);
      fTowerEdges[1] = fEtaBins.at(etaBin);
      fTowerEdges[2] = phiBins.at(phiBin - 1);
      fTowerEdges[3] = phiBins.at(phiBin);

      fECalTowerEnergy = 0.0;
      fHCalTowerEnergy = 0.0;

      fECalTrackEnergy = 0.0;
      fHCalTrackEnergy = 0.0;

      fECalTrackSigma = 0.0;
      fHCalTrackSigma = 0.0;

      fTowerTrackHits = 0;
      fTowerPhotonHits = 0;

      fECalTowerTrackArray->clear();
      fHCalTowerTrackArray->clear();
    }

    // check for track hits
    if(flags & 1)
    {
      ++fTowerTrackHits;

      track = static_cast<Candidate *>(fTrackInputArray->at(number));
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
        fECalTrackEnergy += ecalEnergy;
        ecalSigma = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, momentum.E());
        if(ecalSigma / momentum.E() < track->TrackResolution)
          energyGuess = ecalEnergy;
        else
          energyGuess = momentum.E();

        fECalTrackSigma += (track->TrackResolution) * energyGuess * (track->TrackResolution) * energyGuess;
        fECalTowerTrackArray->emplace_back(track);
      }

      else if(fECalTrackFractions[number] < 1.0E-9 && fHCalTrackFractions[number] > 1.0E-9)
      {
        fHCalTrackEnergy += hcalEnergy;
        hcalSigma = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, momentum.E());
        if(hcalSigma / momentum.E() < track->TrackResolution)
          energyGuess = hcalEnergy;
        else
          energyGuess = momentum.E();

        fHCalTrackSigma += (track->TrackResolution) * energyGuess * (track->TrackResolution) * energyGuess;
        fHCalTowerTrackArray->emplace_back(track);
      }

      else if(fECalTrackFractions[number] < 1.0E-9 && fHCalTrackFractions[number] < 1.0E-9)
      {
        fEFlowTrackOutputArray->emplace_back(track);
      }

      continue;
    }

    // check for photon and electron hits in current tower
    if(flags & 2) ++fTowerPhotonHits;

    particle = static_cast<Candidate *>(fParticleInputArray->at(number));
    momentum = particle->Momentum;
    position = particle->Position;

    // fill current tower
    ecalEnergy = momentum.E() * fECalTowerFractions[number];
    hcalEnergy = momentum.E() * fHCalTowerFractions[number];

    fECalTowerEnergy += ecalEnergy;
    fHCalTowerEnergy += hcalEnergy;

    if(ecalEnergy > fTimingEnergyMin && fTower)
    {
      if(abs(particle->PID) != 11 || !fElectronsFromTrack)
      {
        fTower->ECalEnergyTimePairs.push_back(make_pair<Float_t, Float_t>(ecalEnergy, particle->Position.T()));
      }
    }

    fTower->AddCandidate(particle);
    fTower->Position = position;
  }

  // finalize last tower
  FinalizeTower();
}

//------------------------------------------------------------------------------

void Calorimeter::FinalizeTower()
{
  double energy, pt, eta, phi, r;
  double ecalEnergy, hcalEnergy;
  double ecalNeutralEnergy, hcalNeutralEnergy;

  double ecalSigma, hcalSigma;
  double ecalNeutralSigma, hcalNeutralSigma;

  double weightTrack, weightCalo, bestEnergyEstimate, rescaleFactor;

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

  if(ecalEnergy < fECalEnergyMin || ecalEnergy < fECalEnergySignificanceMin * ecalSigma) ecalEnergy = 0.0;
  if(hcalEnergy < fHCalEnergyMin || hcalEnergy < fHCalEnergySignificanceMin * hcalSigma) hcalEnergy = 0.0;

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

  pt = energy / std::cosh(eta);

  // Time calculation for tower
  fTower->NTimeHits = 0;
  sumWeightedTime = 0.0;
  sumWeight = 0.0;

  for(size_t i = 0; i < fTower->ECalEnergyTimePairs.size(); ++i)
  {
    weight = std::pow((fTower->ECalEnergyTimePairs[i].first), 2);
    sumWeightedTime += weight * fTower->ECalEnergyTimePairs[i].second;
    sumWeight += weight;
    fTower->NTimeHits++;
  }

  // check whether barrel or endcap tower
  if(fTower->Position.Perp() < fTowerRmax && std::fabs(eta) > 0.)
    r = fTower->Position.Z() / std::sinh(eta);
  else
    r = fTower->Position.Pt();

  if(sumWeight > 0.0)
  {
    fTower->Position.SetPtEtaPhiE(r, eta, phi, sumWeightedTime / sumWeight);
  }
  else
  {
    fTower->Position.SetPtEtaPhiE(r, eta, phi, 999999.9);
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
      fPhotonOutputArray->emplace_back(fTower);
    }

    fTowerOutputArray->emplace_back(fTower);
  }

  // fill energy flow candidates
  fECalTrackSigma = std::sqrt(fECalTrackSigma);
  fHCalTrackSigma = std::sqrt(fHCalTrackSigma);

  //compute neutral excesses
  ecalNeutralEnergy = max((ecalEnergy - fECalTrackEnergy), 0.0);
  hcalNeutralEnergy = max((hcalEnergy - fHCalTrackEnergy), 0.0);

  ecalNeutralSigma = ecalNeutralEnergy / std::sqrt(fECalTrackSigma * fECalTrackSigma + ecalSigma * ecalSigma);
  hcalNeutralSigma = hcalNeutralEnergy / std::sqrt(fHCalTrackSigma * fHCalTrackSigma + hcalSigma * hcalSigma);

  // if ecal neutral excess is significant, simply create neutral EflowPhoton tower and clone each track into eflowtrack
  if(ecalNeutralEnergy > fECalEnergyMin && ecalNeutralSigma > fECalEnergySignificanceMin)
  {
    // create new photon tower assuming null mass
    Candidate *tower = static_cast<Candidate *>(fTower->Clone());
    pt = ecalNeutralEnergy / std::cosh(eta);

    tower->Momentum.SetPtEtaPhiE(pt, eta, phi, ecalNeutralEnergy);
    tower->Eem = ecalNeutralEnergy;
    tower->Ehad = 0.0;
    tower->PID = 22;

    fEFlowPhotonOutputArray->emplace_back(tower);

    //clone tracks
    for(Candidate *const &track : *fECalTowerTrackArray)
    {
      Candidate *new_track = static_cast<Candidate *>(track->Clone());
      new_track->AddCandidate(track);

      fEFlowTrackOutputArray->emplace_back(new_track);
    }
  }

  // if neutral excess is not significant, rescale eflow tracks, such that the total charged equals the best measurement given by the calorimeter and tracking
  else if(fECalTrackEnergy > 0.0)
  {
    weightTrack = (fECalTrackSigma > 0.0) ? 1 / (fECalTrackSigma * fECalTrackSigma) : 0.0;
    weightCalo = (ecalSigma > 0.0) ? 1 / (ecalSigma * ecalSigma) : 0.0;

    bestEnergyEstimate = (weightTrack * fECalTrackEnergy + weightCalo * ecalEnergy) / (weightTrack + weightCalo);
    rescaleFactor = bestEnergyEstimate / fECalTrackEnergy;

    //rescale tracks
    for(Candidate *const &track : *fECalTowerTrackArray)
    {
      Candidate *new_track = static_cast<Candidate *>(track->Clone());
      new_track->AddCandidate(track);

      new_track->Momentum *= rescaleFactor;

      fEFlowTrackOutputArray->emplace_back(new_track);
    }
  }

  // if hcal neutral excess is significant, simply create neutral EflowNeutralHadron tower and clone each track into eflowtrack
  if(hcalNeutralEnergy > fHCalEnergyMin && hcalNeutralSigma > fHCalEnergySignificanceMin)
  {
    // create new photon tower
    Candidate *tower = static_cast<Candidate *>(fTower->Clone());
    pt = hcalNeutralEnergy / std::cosh(eta);

    tower->Momentum.SetPtEtaPhiE(pt, eta, phi, hcalNeutralEnergy);
    tower->Ehad = hcalNeutralEnergy;
    tower->Eem = 0.0;

    fEFlowNeutralHadronOutputArray->emplace_back(tower);

    //clone tracks
    for(Candidate *const &track : *fHCalTowerTrackArray)
    {
      Candidate *new_track = static_cast<Candidate *>(track->Clone());
      new_track->AddCandidate(track);

      fEFlowTrackOutputArray->emplace_back(new_track);
    }
  }

  // if neutral excess is not significant, rescale eflow tracks, such that the total charged equals the best measurement given by the calorimeter and tracking
  else if(fHCalTrackEnergy > 0.0)
  {
    weightTrack = (fHCalTrackSigma > 0.0) ? 1 / (fHCalTrackSigma * fHCalTrackSigma) : 0.0;
    weightCalo = (hcalSigma > 0.0) ? 1 / (hcalSigma * hcalSigma) : 0.0;

    bestEnergyEstimate = (weightTrack * fHCalTrackEnergy + weightCalo * hcalEnergy) / (weightTrack + weightCalo);
    rescaleFactor = bestEnergyEstimate / fHCalTrackEnergy;

    //rescale tracks
    for(Candidate *const &track : *fHCalTowerTrackArray)
    {
      Candidate *new_track = static_cast<Candidate *>(track->Clone());
      new_track->AddCandidate(track);
      new_track->Momentum *= rescaleFactor;
      new_track->Momentum.SetPtEtaPhiM(track->Momentum.Pt() * rescaleFactor, track->Momentum.Eta(), track->Momentum.Phi(), track->Momentum.M());

      fEFlowTrackOutputArray->emplace_back(new_track);
    }
  }
}

//------------------------------------------------------------------------------

double Calorimeter::LogNormal(Double_t mean, Double_t sigma)
{
  if(mean > 0.0)
  {
    const double b = std::sqrt(std::log((1.0 + (sigma * sigma) / (mean * mean)))),
                 a = std::log(mean) - 0.5 * b * b;

    return std::exp(a + b * gRandom->Gaus(0.0, 1.0));
  }
  else
    return 0.0;
}

//------------------------------------------------------------------------------

REGISTER_MODULE("Calorimeter", Calorimeter);
