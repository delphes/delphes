
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
 *  \author A. Chattopadhyay - UPRM
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

#include <set>

using namespace std;

class SimpleCalorimeter: public DelphesModule
{
public:
  SimpleCalorimeter(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    //
    fEnergyMin(Steer<double>("EnergyMin", 0.0)), // read min E value for towers to be saved
    fEnergySignificanceMin(Steer<double>("EnergySignificanceMin", 0.0)),
    fIsEcal(Steer<bool>("IsEcal", false)), // flag that says if current calo is Ecal of Hcal (will then fill correct values of Eem and Ehad)
    fSmearTowerCenter(Steer<bool>("SmearTowerCenter", true)), // switch on or off the dithering of the center of calorimeter towers
    fPhiBins(Steer<std::unordered_map<double, std::vector<double> > >("EtaPhiBins")),
    fFractionMap(Steer<TFractionMap>("EnergyFraction")), // energy fractions for different particles
    //
    fResolutionFormula(std::make_unique<DelphesFormula>())
  {
    for(const std::pair<double, std::vector<double> > etaPhiBins : fPhiBins) // auto would avoid a copy
      fEtaBins.emplace_back(etaPhiBins.first);
    std::sort(fEtaBins.begin(), fEtaBins.end());

    // for blind calorimeter: read insensitive bins (convert centers -> indices)
    // Loop over blind eta-phi pairs
    for(const auto &[etaEdge, phiEdge] :
      Steer<std::vector<std::pair<double, double> > >("InsensitiveEtaPhiBins"))
    {
      // Find closest eta bin index (using bin centers)
      auto itEta = std::min_element(fEtaBins.begin(), fEtaBins.end(), [etaEdge](double a, double b) { return std::abs(a - etaEdge) < std::abs(b - etaEdge); });
      if(itEta == fEtaBins.end()) continue;
      const short etaBinIndex = static_cast<Short_t>(std::distance(fEtaBins.begin(), itEta));

      // Find closest phi bin index (using bin centers)
      const std::vector<double> &phiVec = fPhiBins.at(fEtaBins.at(etaBinIndex));
      if(phiVec.empty()) continue;

      auto itPhi = std::min_element(phiVec.begin(), phiVec.end(), [phiEdge](double a, double b) { return std::abs(a - phiEdge) < std::abs(b - phiEdge); });
      if(itPhi == phiVec.end()) continue;
      const short phiBinIndex = static_cast<Short_t>(std::distance(phiVec.begin(), itPhi));

      // Insert bin into insensitive set
      fInsensitiveBinSet.insert(std::make_pair(etaBinIndex, phiBinIndex));
    }

    // set default energy fraction value
    fFractionMap[0] = 1.;

    // read resolution formulas
    fResolutionFormula->Compile(Steer<std::string>("ResolutionFormula", "0"));
  }

  void Init() override
  {
    fParticleInputArray = ImportArray(Steer<std::string>("ParticleInputArray", "ParticlePropagator/particles"));
    fTrackInputArray = ImportArray(Steer<std::string>("TrackInputArray", "ParticlePropagator/tracks"));
    fTowerOutputArray = ExportArray(Steer<std::string>("TowerOutputArray", "towers"));
    fEFlowTrackOutputArray = ExportArray(Steer<std::string>("EFlowTrackOutputArray", "eflowTracks"));
    fEFlowTowerOutputArray = ExportArray(Steer<std::string>("EFlowTowerOutputArray", "eflowTowers"));
  }
  void Process() override;

private:
  typedef std::map<long long, double> TFractionMap; //!
  typedef std::map<double, std::set<double> > TBinMap; //!

  void FinalizeTower();
  double LogNormal(double mean, double sigma);

  const double fEnergyMin;
  const double fEnergySignificanceMin;

  const bool fIsEcal; //!
  const bool fSmearTowerCenter;

  const std::unordered_map<double, std::vector<double> > fPhiBins;
  std::vector<double> fEtaBins;

  Candidate *fTower{nullptr};
  double fTowerEta, fTowerPhi, fTowerEdges[4];

  double fTowerEnergy;
  double fNeutralEnergy;
  double fTrackEnergy;

  double fTowerEnergyFromPU;
  double fTrackEnergyFromPU;
  double fNeutralEnergyFromPU;

  double fTowerTime;
  double fTrackTime;

  double fTowerRmax;

  double fTowerTimeWeight;
  double fTrackTimeWeight;

  int fTowerTrackHits, fTowerPhotonHits;

  double fTrackSigma;

  TFractionMap fFractionMap; //!
  TBinMap fBinMap; //!

  std::vector<unsigned long long> fTowerHits;

  std::vector<double> fTowerFractions;
  std::vector<double> fTrackFractions;

  // Insensitive bins stored as integer indices (etaBin, phiBin)
  std::set<std::pair<Short_t, Short_t> > fInsensitiveBinSet;

  // Flag telling whether the *current* tower is insensitive or not
  inline bool IsTowerInsensitive(Short_t etaBin, Short_t phiBin) const
  {
    return fInsensitiveBinSet.find(std::make_pair(etaBin, phiBin)) != fInsensitiveBinSet.end();
  }

  const std::unique_ptr<DelphesFormula> fResolutionFormula; //!

  CandidatesCollection fParticleInputArray; //!
  CandidatesCollection fTrackInputArray; //!

  CandidatesCollection fTowerOutputArray; //!

  CandidatesCollection fEFlowTrackOutputArray; //!
  CandidatesCollection fEFlowTowerOutputArray; //!

  CandidatesCollection fTowerTrackArray; //!
};

//------------------------------------------------------------------------------

void SimpleCalorimeter::Process()
{
  fTowerOutputArray->clear();
  fEFlowTrackOutputArray->clear();
  fEFlowTowerOutputArray->clear();

  TLorentzVector position, momentum;
  Short_t phiBin, flags;
  int number;
  unsigned long long towerHit, towerEtaPhi, hitEtaPhi;
  double fraction;
  double energy;
  double sigma;
  double energyGuess;

  int pdgCode;

  TFractionMap::iterator itFractionMap;

  vector<unsigned long long>::iterator itTowerHits;

  DelphesFactory *factory = GetFactory();
  fTowerHits.clear();
  fTowerFractions.clear();
  fTrackFractions.clear();

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
      itFractionMap = fFractionMap.find(0);

    fraction = itFractionMap->second;
    fTowerFractions.push_back(fraction);

    if(fraction < 1.0E-9) continue;

    // find eta bin [1, fEtaBins.size - 1]
    std::vector<double>::iterator itEtaBin = lower_bound(fEtaBins.begin(), fEtaBins.end(), particlePosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    const short etaBin = std::distance(fEtaBins.begin(), itEtaBin);

    // phi bins for given eta bin
    const std::vector<double> &phiBins = fPhiBins.at(*itEtaBin);

    // find phi bin [1, phiBins.size - 1]
    std::vector<double>::const_iterator itPhiBin = std::lower_bound(phiBins.begin(), phiBins.end(), particlePosition.Phi());
    if(itPhiBin == phiBins.begin() || itPhiBin == phiBins.end()) continue;
    phiBin = std::distance(phiBins.begin(), itPhiBin);

    flags = 0;
    flags |= (pdgCode == 11 || pdgCode == 22) << 1;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for particle number}
    towerHit = ((unsigned long long)(etaBin) << 48) | ((unsigned long long)(phiBin) << 32) | ((unsigned long long)(flags) << 24) | (unsigned long long)(number);

    fTowerHits.push_back(towerHit);

    // skip insensitive calorimeter bins entirely for particles
    if(IsTowerInsensitive(etaBin, phiBin))
    {
      fTower = nullptr;
    }
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

    fraction = itFractionMap->second;

    fTrackFractions.push_back(fraction);

    // find eta bin [1, fEtaBins.size - 1]
    std::vector<double>::iterator itEtaBin = std::lower_bound(fEtaBins.begin(), fEtaBins.end(), trackPosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    const short etaBin = std::distance(fEtaBins.begin(), itEtaBin);

    // phi bins for given eta bin
    const std::vector<double> &phiBins = fPhiBins.at(*itEtaBin);

    // find phi bin [1, phiBins.size - 1]
    std::vector<double>::const_iterator itPhiBin = std::lower_bound(phiBins.begin(), phiBins.end(), trackPosition.Phi());
    if(itPhiBin == phiBins.begin() || itPhiBin == phiBins.end()) continue;
    const short phiBin = std::distance(phiBins.cbegin(), itPhiBin);

    flags = 1;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for track number}
    towerHit = ((unsigned long long)(etaBin) << 48) | ((unsigned long long)(phiBin) << 32) | ((unsigned long long)(flags) << 24) | (unsigned long long)(number);

    fTowerHits.push_back(towerHit);
    // skip insensitive calo bins entirely for particles
    if(IsTowerInsensitive(etaBin, phiBin))
    {
      fTower = nullptr;
    }
  }

  // all hits are sorted first by eta bin number, then by phi bin number,
  // then by flags and then by particle or track number
  sort(fTowerHits.begin(), fTowerHits.end());

  // loop over all hits
  towerEtaPhi = 0;
  fTower = nullptr;
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
      if(fTower) FinalizeTower();

      // create new tower
      fTower = factory->NewCandidate();
      const short phiBin = (towerHit >> 32) & 0x000000000000FFFFLL;
      const short etaBin = (towerHit >> 48) & 0x000000000000FFFFLL;

      //mark fTower nullptr
      if(IsTowerInsensitive(etaBin, phiBin))
      {
        fTower = nullptr; // insensitive tower: no creation
        // do NOT continue here! preserve hit loop for ordering
      }

      // phi bins for given eta bin
      const std::vector<double> &phiBins = fPhiBins.at(fEtaBins.at(etaBin));

      // calculate eta and phi of the tower's center
      fTowerEta = 0.5 * (fEtaBins.at(etaBin - 1) + fEtaBins.at(etaBin));
      fTowerPhi = 0.5 * (phiBins.at(phiBin - 1) + phiBins.at(phiBin));

      fTowerEdges[0] = fEtaBins.at(etaBin - 1);
      fTowerEdges[1] = fEtaBins.at(etaBin);
      fTowerEdges[2] = phiBins.at(phiBin - 1);
      fTowerEdges[3] = phiBins.at(phiBin);

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

      fTowerTrackArray->clear();
    }

    // We already handled it above; fTower == nullptr if insensitive
    // This prevents skipping hits that are already stored in fTowerHits

    // check for track hits
    if(flags & 1)
    {
      ++fTowerTrackHits;
      Candidate *track = static_cast<Candidate *>(fTrackInputArray->at(number));
      momentum = track->Momentum;
      position = track->Position;

      energy = momentum.E() * fTrackFractions[number];

      fTrackTime += std::sqrt(energy) * position.T();
      fTrackTimeWeight += std::sqrt(energy);

      if(fTrackFractions[number] > 1.0E-9)
      {

        // compute total charged energy
        fTrackEnergy += energy;

        // compute energy contribution from PU tracks
        if(track->IsPU) fTrackEnergyFromPU += energy;

        sigma = fResolutionFormula->Eval(0.0, fTowerEta, 0.0, momentum.E());
        if(sigma / momentum.E() < track->TrackResolution)
          energyGuess = energy;
        else
          energyGuess = momentum.E();

        fTrackSigma += ((track->TrackResolution) * energyGuess) * ((track->TrackResolution) * energyGuess);
        if(fTower) fTowerTrackArray->emplace_back(track); // add only if tower exists
      }
      else
      {
        if(fTower) fEFlowTrackOutputArray->emplace_back(track);
      }
      continue;
    }

    // check for photon and electron hits in current tower
    if(flags & 2) ++fTowerPhotonHits;

    Candidate *particle = static_cast<Candidate *>(fParticleInputArray->at(number));
    momentum = particle->Momentum;
    position = particle->Position;

    // fill current tower
    energy = momentum.E() * fTowerFractions[number];

    if(!(flags & 1))
    {
      // compute total neutral energy
      fNeutralEnergy += energy;
      if(particle->IsPU) fNeutralEnergyFromPU += energy;
    }

    if(fTower) // add only if tower exists
    {
      fTowerEnergy += energy;
      fTowerTime += energy * energy * position.T(); //sigma_t ~ 1/E
      fTowerTimeWeight += energy * energy;
      fTower->AddCandidate(particle);
      fTower->Position = position;
      if(particle->IsPU) fTowerEnergyFromPU += energy;
    }
  }

  // finalize last tower (only if it was sensitive)
  if(fTower) FinalizeTower();
}

//------------------------------------------------------------------------------

void SimpleCalorimeter::FinalizeTower()
{
  double energy, neutralEnergy, pt, eta, phi, r;
  double sigma, neutralSigma;
  double time;

  double weightTrack, weightCalo, bestEnergyEstimate, rescaleFactor;

  double hardEnergyFraction;

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

  pt = energy / std::cosh(eta);

  // endcap
  if(std::fabs(fTower->Position.Pt() - fTowerRmax) > 1.e-06 && std::fabs(eta) > 0.)
  {
    r = fTower->Position.Z() / std::sinh(eta);
  }
  // barrel
  else
  {
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
  if(energy > 0.0) fTowerOutputArray->emplace_back(fTower);

  // e-flow candidates

  //compute neutral excess

  fTrackSigma = std::sqrt(fTrackSigma);
  neutralEnergy = max((energy - fTrackEnergy), 0.0);

  //compute sigma_trk total
  neutralSigma = neutralEnergy / std::sqrt(fTrackSigma * fTrackSigma + sigma * sigma);

  // if neutral excess is significant, simply create neutral Eflow tower and clone each track into eflowtrack
  if(neutralEnergy > fEnergyMin && neutralSigma > fEnergySignificanceMin)
  {
    // create new photon tower
    Candidate *tower = static_cast<Candidate *>(fTower->Clone());
    pt = neutralEnergy / std::cosh(eta);

    tower->Eem = (!fIsEcal) ? 0 : neutralEnergy;
    tower->Ehad = (fIsEcal) ? 0 : neutralEnergy;
    tower->PID = (fIsEcal) ? 22 : 0;

    tower->Momentum.SetPtEtaPhiE(pt, eta, phi, neutralEnergy);

    // compute hard energy fraction for tower
    hardEnergyFraction = (fNeutralEnergy > 1.0E-09) ? (fNeutralEnergy - fNeutralEnergyFromPU) / (fNeutralEnergy) : 1.0;

    if(hardEnergyFraction < 0.5)
      tower->IsPU = 1;
    else
      tower->IsPU = 0;

    tower->BetaStar = hardEnergyFraction;

    fEFlowTowerOutputArray->emplace_back(tower);

    for(Candidate *const &track : *fTowerTrackArray)
    {
      Candidate *new_track = static_cast<Candidate *>(track->Clone());
      new_track->AddCandidate(track);

      fEFlowTrackOutputArray->emplace_back(new_track);
    }
  }

  // if neutral excess is not significant, rescale eflow tracks, such that the total charged equals the best measurement given by the calorimeter and tracking
  else if(fTrackEnergy > 0.0)
  {

    weightTrack = (fTrackSigma > 0.0) ? 1 / (fTrackSigma * fTrackSigma) : 0.0;
    weightCalo = (sigma > 0.0) ? 1 / (sigma * sigma) : 0.0;

    bestEnergyEstimate = (weightTrack * fTrackEnergy + weightCalo * energy) / (weightTrack + weightCalo);
    rescaleFactor = bestEnergyEstimate / fTrackEnergy;

    for(Candidate *const &track : *fTowerTrackArray)
    {
      Candidate *new_track = static_cast<Candidate *>(track->Clone());
      new_track->AddCandidate(track);
      new_track->Momentum.SetPtEtaPhiM(new_track->Momentum.Pt() * rescaleFactor, new_track->Momentum.Eta(), new_track->Momentum.Phi(), new_track->Momentum.M());
      fEFlowTrackOutputArray->emplace_back(new_track);
    }
  }
}

//------------------------------------------------------------------------------

double SimpleCalorimeter::LogNormal(double mean, double sigma)
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

REGISTER_MODULE("SimpleCalorimeter", SimpleCalorimeter);
