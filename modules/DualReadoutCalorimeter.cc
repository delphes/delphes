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

/** \class DualReadoutCalorimeter
 *
 *  Fills DualReadoutCalorimeter towers, performs DualReadoutCalorimeter resolution smearing,
 *  and creates energy flow objects (tracks, photons, and neutral hadrons).
 *
 *  \author M. Selvaggi - CERN
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

#include <iostream>
#include <map>
#include <set>
#include <vector>

using namespace std;

class DualReadoutCalorimeter: public DelphesModule
{
public:
  explicit DualReadoutCalorimeter(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fTimingEnergyMin(Steer<double>("TimingEnergyMin", 4.)), // read min E value for timing measurement in ECAL
    // should be optimised depending calorimeter resolution
    fECalMinSignificance(Steer<double>("ECalMinSignificance", 0.0)),
    fHCalMinSignificance(Steer<double>("HCalMinSignificance", 0.0)),
    // switch on or off the dithering of the center of DualReadoutCalorimeter towers
    fSmearTowerCenter(Steer<bool>("SmearTowerCenter", true)),
    // switch on or off the log normal smearing (gaussian if false)
    fSmearLogNormal(Steer<bool>("SmearLogNormal", true)),
    fPhiBins(Steer<std::unordered_map<double, std::vector<double> > >("EtaPhiBins")),
    fFractionMap(Steer<TFractionMap>("EnergyFraction")), // ECAL/HCAL energy fractions for different particles
    fECalResolutionFormula(std::make_unique<DelphesFormula>()),
    fHCalResolutionFormula(std::make_unique<DelphesFormula>()),
    fTowerTrackArray(std::make_shared<std::vector<Candidate *> >())
  {
    for(const std::pair<double, std::vector<double> > etaPhiBins : fPhiBins) // auto would avoid a copy
      fEtaBins.emplace_back(etaPhiBins.first);
    std::sort(fEtaBins.begin(), fEtaBins.end());

    // set default energy fractions values
    fFractionMap[0] = std::make_pair(0.0, 1.0);

    // read resolution formulas
    fECalResolutionFormula->Compile(Steer<std::string>("ECalResolutionFormula", "0"));
    fHCalResolutionFormula->Compile(Steer<std::string>("HCalResolutionFormula", "0"));
  }

  void Init() override
  {
    // import arrays with output from other modules
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
  typedef std::map<unsigned long long, std::pair<double, double> > TFractionMap; //!
  typedef std::map<double, std::set<double> > TBinMap; //!

  // For timing
  // So far this flag needs to be false
  // Curved extrapolation not supported
  const bool fElectronsFromTrack{false};
  const double fTimingEnergyMin;

  const double fECalMinSignificance;
  const double fHCalMinSignificance;

  const bool fSmearTowerCenter;
  const bool fSmearLogNormal;

  const std::unordered_map<double, std::vector<double> > fPhiBins;
  std::vector<double> fEtaBins;

  TFractionMap fFractionMap; //!

  const std::unique_ptr<DelphesFormula> fECalResolutionFormula; //!
  const std::unique_ptr<DelphesFormula> fHCalResolutionFormula; //!

  Candidate *fTower{nullptr};
  double fTowerEta, fTowerPhi, fTowerEdges[4];
  double fECalTowerEnergy, fHCalTowerEnergy;
  double fECalTrackEnergy, fHCalTrackEnergy;
  double fTrackEnergy;
  double fTowerRmax;

  int fTowerTrackHits, fTowerPhotonHits;

  double fTrackSigma;

  double fTowerTime;
  double fTowerTimeWeight;

  std::vector<unsigned long long> fTowerHits;

  std::vector<double> fECalTowerFractions;
  std::vector<double> fHCalTowerFractions;

  std::vector<double> fECalTrackFractions;
  std::vector<double> fHCalTrackFractions;

  void FinalizeTower();
  double LogNormal(double mean, double sigma);
  double TruncatedGaussian(double mean, double sigma);

  //CandidatesCollection fECalTowerTrackArray; //!
  //CandidatesCollection fHCalTowerTrackArray; //!
  const CandidatesCollection fTowerTrackArray; //!

  CandidatesCollection fParticleInputArray; //!
  CandidatesCollection fTrackInputArray; //!

  CandidatesCollection fTowerOutputArray; //!
  CandidatesCollection fPhotonOutputArray; //!

  CandidatesCollection fEFlowTrackOutputArray; //!
  CandidatesCollection fEFlowPhotonOutputArray; //!
  CandidatesCollection fEFlowNeutralHadronOutputArray; //!
};

//------------------------------------------------------------------------------

void DualReadoutCalorimeter::Process()
{
  fTowerOutputArray->clear();
  fPhotonOutputArray->clear();
  fEFlowTrackOutputArray->clear();
  fEFlowPhotonOutputArray->clear();
  fEFlowNeutralHadronOutputArray->clear();

  DelphesFactory *factory = GetFactory();
  fTowerHits.clear();
  fECalTowerFractions.clear();
  fHCalTowerFractions.clear();
  fECalTrackFractions.clear();
  fHCalTrackFractions.clear();

  // loop over all particles
  size_t number = 0;
  fTowerRmax = 0.;

  //cout<<"--------- new event ---------- "<<endl;
  for(Candidate *const &particle : *fParticleInputArray)
  {
    const TLorentzVector &particlePosition = particle->Position;

    // compute maximum radius (needed in FinalizeTower to assess whether barrel or endcap tower)
    if(particlePosition.Perp() > fTowerRmax)
      fTowerRmax = particlePosition.Perp();

    const int pdgCode = std::abs(particle->PID);

    TFractionMap::iterator itFractionMap = fFractionMap.find(pdgCode);
    if(itFractionMap == fFractionMap.end())
      itFractionMap = fFractionMap.find(0);

    const double ecalFraction = itFractionMap->second.first;
    const double hcalFraction = itFractionMap->second.second;

    fECalTowerFractions.push_back(ecalFraction);
    fHCalTowerFractions.push_back(hcalFraction);

    if(ecalFraction < 1.0E-9 && hcalFraction < 1.0E-9) continue;

    // find eta bin [1, fEtaBins.size - 1]
    std::vector<double>::const_iterator itEtaBin = std::lower_bound(fEtaBins.begin(), fEtaBins.end(), particlePosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    const short etaBin = std::distance(fEtaBins.cbegin(), itEtaBin);

    // phi bins for given eta bin
    const std::vector<double> &phiBins = fPhiBins.at(*itEtaBin);

    // find phi bin [1, phiBins.size - 1]
    std::vector<double>::const_iterator itPhiBin = std::lower_bound(phiBins.begin(), phiBins.end(), particlePosition.Phi());
    if(itPhiBin == phiBins.begin() || itPhiBin == phiBins.end()) continue;
    const short phiBin = std::distance(phiBins.cbegin(), itPhiBin);

    short flags = 0;
    flags |= (pdgCode == 11 || pdgCode == 22) << 1;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for particle number}
    unsigned long long towerHit = ((unsigned long long)(etaBin) << 48) | ((unsigned long long)(phiBin) << 32) | ((unsigned long long)(flags) << 24) | (unsigned long long)(number);

    fTowerHits.push_back(towerHit);
    ++number;
  }

  // loop over all tracks
  number = 0;
  for(Candidate *const &track : *fTrackInputArray)
  {
    const TLorentzVector &trackPosition = track->Position;

    const int pdgCode = std::abs(track->PID);

    TFractionMap::iterator itFractionMap = fFractionMap.find(pdgCode);
    if(itFractionMap == fFractionMap.end())
      itFractionMap = fFractionMap.find(0);

    const double ecalFraction = itFractionMap->second.first;
    const double hcalFraction = itFractionMap->second.second;

    fECalTrackFractions.push_back(ecalFraction);
    fHCalTrackFractions.push_back(hcalFraction);

    // find eta bin [1, fEtaBins.size - 1]
    std::vector<double>::const_iterator itEtaBin = std::lower_bound(fEtaBins.begin(), fEtaBins.end(), trackPosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    const short etaBin = std::distance(fEtaBins.cbegin(), itEtaBin);

    // phi bins for given eta bin
    const std::vector<double> &phiBins = fPhiBins.at(*itEtaBin);

    // find phi bin [1, phiBins.size - 1]
    std::vector<double>::const_iterator itPhiBin = std::lower_bound(phiBins.begin(), phiBins.end(), trackPosition.Phi());
    if(itPhiBin == phiBins.begin() || itPhiBin == phiBins.end()) continue;
    const short phiBin = std::distance(phiBins.cbegin(), itPhiBin);

    const short flags = 1;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for track number}
    unsigned long long towerHit = ((unsigned long long)(etaBin) << 48) | ((unsigned long long)(phiBin) << 32) | ((unsigned long long)(flags) << 24) | (unsigned long long)(number);

    fTowerHits.push_back(towerHit);
    ++number;
  }

  // all hits are sorted first by eta bin number, then by phi bin number,
  // then by flags and then by particle or track number
  std::sort(fTowerHits.begin(), fTowerHits.end());

  // loop over all hits
  unsigned long long towerEtaPhi = 0;
  fTower = 0;
  for(const unsigned long long &towerHit : fTowerHits)
  {
    const short flags = (towerHit >> 24) & 0x00000000000000FFLL;
    const short number = (towerHit) & 0x0000000000FFFFFFLL;
    const unsigned long long hitEtaPhi = towerHit >> 32;

    if(towerEtaPhi != hitEtaPhi)
    {
      // switch to next tower
      towerEtaPhi = hitEtaPhi;

      // finalize previous tower
      FinalizeTower();

      // create new tower
      fTower = factory->NewCandidate();

      const short phiBin = (towerHit >> 32) & 0x000000000000FFFFLL;
      const short etaBin = (towerHit >> 48) & 0x000000000000FFFFLL;

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
      fTrackEnergy = 0.0;

      fTrackSigma = 0.0;

      fTowerTrackHits = 0;
      fTowerPhotonHits = 0;

      fTowerTime = 0.0;
      fTowerTimeWeight = 0.0;

      //fECalTowerTrackArray->clear();
      //fHCalTowerTrackArray->clear();
      fTowerTrackArray->clear();
    }

    // check for track hits
    if(flags & 1)
    {
      ++fTowerTrackHits;

      Candidate *track = static_cast<Candidate *>(fTrackInputArray->at(number));
      const TLorentzVector &momentum = track->Momentum;

      const double ecalEnergy = momentum.E() * fECalTrackFractions[number];
      const double hcalEnergy = momentum.E() * fHCalTrackFractions[number];
      const double energy = ecalEnergy + hcalEnergy;

      if(ecalEnergy > fTimingEnergyMin && fTower)
      {
        if(fElectronsFromTrack)
        {
          fTower->ECalEnergyTimePairs.push_back(make_pair<Float_t, Float_t>(ecalEnergy, track->Position.T()));
        }
      }

      // in Dual Readout we do not care if tracks are ECAL of HCAL
      if(fECalTrackFractions[number] > 1.0E-9 || fHCalTrackFractions[number] > 1.0E-9)
      {
        fTrackEnergy += energy;
        // this sigma will be used to determine whether neutral excess is significant. We choose the resolution according to bthe higest deposited fraction (in practice had for charged hadrons and em for electrons)
        double sigma = 0.;
        if(fHCalTrackFractions[number] > 0)
          sigma = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, momentum.E());
        else
          sigma = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, momentum.E());

        double energyGuess = 0.;
        if(sigma / momentum.E() < track->TrackResolution)
          energyGuess = ecalEnergy + hcalEnergy;
        else
          energyGuess = momentum.E();

        fTrackSigma += (track->TrackResolution) * energyGuess * (track->TrackResolution) * energyGuess;
        fTowerTrackArray->emplace_back(track);
      }
      else
      {
        fEFlowTrackOutputArray->emplace_back(track);
      }

      continue;
    }

    // check for photon and electron hits in current tower
    if(flags & 2) ++fTowerPhotonHits;

    Candidate *particle = static_cast<Candidate *>(fParticleInputArray->at(number));
    const TLorentzVector &momentum = particle->Momentum;
    const TLorentzVector &position = particle->Position;

    // fill current tower
    const double ecalEnergy = momentum.E() * fECalTowerFractions[number];
    const double hcalEnergy = momentum.E() * fHCalTowerFractions[number];

    fECalTowerEnergy += ecalEnergy;
    fHCalTowerEnergy += hcalEnergy;

    // assume combined timing measurements in ECAL/HCAL sections
    fTowerTime += (ecalEnergy + hcalEnergy) * position.T(); //sigma_t ~ 1/sqrt(E)
    fTowerTimeWeight += ecalEnergy + hcalEnergy;

    fTower->AddCandidate(particle);
    fTower->Position = position;
  }

  // finalize last tower
  FinalizeTower();
}

//------------------------------------------------------------------------------

void DualReadoutCalorimeter::FinalizeTower()
{
  fTowerOutputArray->clear();
  fPhotonOutputArray->clear();
  fEFlowTrackOutputArray->clear();
  fEFlowPhotonOutputArray->clear();
  fEFlowNeutralHadronOutputArray->clear();

  double eta, phi, r, time;
  double neutralEnergy;

  double caloSigma, trackCaloSigma;
  double neutralSignificance;
  double neutralMinPFSignificance;

  double weightTrack, weightCalo, bestEnergyEstimate, rescaleFactor;
  Bool_t isPureEM = false;

  Bool_t debug = false;
  if(!fTower) return;
  if(debug) cout << "-----------------------------------------------------------------------" << endl;
  if(debug) cout << "New Tower: " << fECalTowerEnergy << "," << fHCalTowerEnergy << "," << fHCalTowerEnergy << "," << fTowerEta << endl;

  if(debug) cout << "   gen particles in tower :" << fTower->GetCandidates().size() << endl;
  for(Candidate *const &candidate : fTower->GetCandidates())
  {
    //cout<<": " << <<endl;
    TLorentzVector mom = candidate->Momentum;
    if(debug) cout << "      gen particle: " << candidate->PID << "," << mom.E() << "," << mom.Eta() << "," << mom.Phi() << endl;
  }

  double energy = 0.;
  // if no hadronic energy, use ECAL resolution
  if(fHCalTowerEnergy <= 0)
  {
    energy = fECalTowerEnergy;
    caloSigma = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, energy);
    isPureEM = true;
    if(debug) cout << "   using ECAL energy: " << energy << ", " << caloSigma << endl;
  }

  // if hadronic fraction > 0, use HCAL resolution
  else
  {
    energy = fECalTowerEnergy + fHCalTowerEnergy;
    caloSigma = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, energy);
    if(debug) cout << "   using HCAL energy: " << energy << ", " << caloSigma << endl;
  }

  if(fSmearLogNormal)
    energy = LogNormal(energy, caloSigma);
  else
    //energy = TruncatedGaussian(energy, caloSigma);
    energy = gRandom->Gaus(energy, caloSigma);

  if(debug) cout << "   smeared energy: " << energy << endl;

  if(energy < 0.) energy = 0.;
  // set tower energy to 0 when energy deposit is not significant

  if(isPureEM)
  {
    // estimate resolution from the measurement this time
    caloSigma = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, energy);
    energy = (energy > fECalMinSignificance * caloSigma) ? energy : 0.;
  }
  else
  {
    // estimate resolution from the measurement this time
    caloSigma = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, energy);
    energy = (energy > fHCalMinSignificance * caloSigma) ? energy : 0.;
  }

  if(debug) cout << "   smeared energy: " << energy << endl;

  // ---------------------------------------------------------------------------
  // compute calo tower properties
  // ---------------------------------------------------------------------------

  time = (fTowerTimeWeight < 1.0E-09) ? 0.0 : fTowerTime / fTowerTimeWeight;

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

  // check whether barrel or endcap tower

  if(std::fabs(fTower->Position.Pt() - fTowerRmax) > 1.e-06 && std::fabs(eta) > 0.) // endcap
    r = fTower->Position.Z() / std::sinh(eta);
  else // barrel
    r = fTower->Position.Pt();

  time = (fTowerTimeWeight < 1.0E-09) ? 0.0 : fTowerTime / fTowerTimeWeight;

  fTower->Position.SetPtEtaPhiE(r, eta, phi, time);
  fTower->L = fTower->Position.Vect().Mag();
  fTower->Etrk = fTrackEnergy;

  // these are stored for debugging purposes, should not be used since they are
  // based on MC truth
  fTower->Eem = fECalTowerEnergy;
  fTower->Ehad = fHCalTowerEnergy;

  fTower->Edges[0] = fTowerEdges[0];
  fTower->Edges[1] = fTowerEdges[1];
  fTower->Edges[2] = fTowerEdges[2];
  fTower->Edges[3] = fTowerEdges[3];

  if(isPureEM)
  {
    // assume massless photon hypothesis
    fTower->PID = 22;
    const double pt = energy / std::cosh(eta);
    fTower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
  }
  // if hadronic fraction > 0, use HCAL resolution
  else
  {

    // assume pion hypothesis for hadronic deposit. This can be corrected later by accessing particle energy in the output
    fTower->PID = 211;
    const double mass = 0.13957;
    const double p = (energy > mass) ? std::sqrt(energy * energy - mass * mass) : 0.;
    const double pt = p / std::cosh(eta);
    fTower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
  }

  if(energy > 0.0)
  {
    if(fTowerPhotonHits > 0 && fTowerTrackHits == 0)
    {
      fPhotonOutputArray->emplace_back(fTower);
    }

    if(debug) cout << "   creating tower with energy: " << energy << endl;
    if(debug) cout << "   creating tower with PID: " << fTower->PID << endl;
    if(debug) cout << "   creating tower with track energy: " << fTower->Etrk << endl;

    fTowerOutputArray->emplace_back(fTower);
  }

  // ---------------------------------------------------------------------------
  // now do particle-flow
  // ---------------------------------------------------------------------------

  fTrackSigma = std::sqrt(fTrackSigma);
  neutralEnergy = max((energy - fTrackEnergy), 0.0);

  if(isPureEM)
  {
    neutralMinPFSignificance = fECalMinSignificance;
  }
  else
  {
    neutralMinPFSignificance = fHCalMinSignificance;
  }

  // combined track calo resolution
  trackCaloSigma = std::sqrt(fTrackSigma * fTrackSigma + caloSigma * caloSigma);
  neutralSignificance = (trackCaloSigma > 0) ? neutralEnergy / trackCaloSigma : 0.;

  if(debug) cout << "Doing PF here: " << endl;
  if(debug) cout << "   track energy: " << fTrackEnergy << endl;
  if(debug) cout << "   calo energy: " << energy << endl;
  if(debug) cout << "   neutral energy: " << neutralEnergy << endl;
  if(debug) cout << "   track sigma: " << fTrackSigma << endl;
  if(debug) cout << "   calo sigma: " << caloSigma << endl;
  if(debug) cout << "   track calo sigma: " << trackCaloSigma << endl;
  if(debug) cout << "   neutral significance: " << neutralSignificance << endl;

  // now do case where at least one track points to tower and the nuetral excess is signficant
  // i.e pi+ and neutron or electron and photon hitting same tower
  if(neutralSignificance > neutralMinPFSignificance)
  {

    Candidate *tower = static_cast<Candidate *>(fTower->Clone());
    if(isPureEM)
    {
      tower->Eem = neutralEnergy;
      tower->Ehad = 0.0;
      tower->PID = 22;
      const double pt = neutralEnergy / std::cosh(eta);
      tower->Momentum.SetPtEtaPhiE(pt, eta, phi, neutralEnergy);
      fEFlowPhotonOutputArray->emplace_back(tower);
    }
    else
    {
      tower->Eem = 0;
      tower->Ehad = neutralEnergy;
      tower->PID = 130;
      const double mass = 0.497611;
      const double p = (neutralEnergy > mass) ? std::sqrt(neutralEnergy * neutralEnergy - mass * mass) : 0.;
      const double pt = p / std::cosh(eta);
      if(p > 0)
      {
        tower->Momentum.SetPtEtaPhiE(pt, eta, phi, neutralEnergy);
        fEFlowNeutralHadronOutputArray->emplace_back(tower);
      }
    }

    if(debug) cout << "       creating neutral excess with energy, eta, phi: " << neutralEnergy << "," << eta << "," << phi << endl;
    if(debug) cout << "       creating neutral excess with PID: " << tower->PID << endl;
    if(debug) cout << "       creating neutral excess with track energy: " << tower->Etrk << endl;
    if(debug) cout << "       " << endl;

    // now clone tracks
    for(Candidate *const &track : *fTowerTrackArray)
    {
      Candidate *new_track = static_cast<Candidate *>(track->Clone());
      new_track->AddCandidate(track);
      fEFlowTrackOutputArray->emplace_back(new_track);
    }
  }

  // now do case where at track points to tower and the neutral excess is NOT signficant
  // if neutral excess is not significant, rescale eflow tracks, such that the total
  // charged equals the best measurement given by the DualReadoutCalorimeter and tracking
  else if(fTrackEnergy > 0)
  {
    if(debug) cout << "   no significant neutral excess found:" << endl;
    if(debug) cout << "   neutral energy: " << energy << ", " << fTrackEnergy << ", " << neutralEnergy << endl;
    if(debug) cout << "       " << endl;

    weightTrack = (fTrackSigma > 0.0) ? 1 / (fTrackSigma * fTrackSigma) : 0.0;
    weightCalo = (caloSigma > 0.0) ? 1 / (caloSigma * caloSigma) : 0.0;

    bestEnergyEstimate = (weightTrack * fTrackEnergy + weightCalo * energy) / (weightTrack + weightCalo);
    rescaleFactor = bestEnergyEstimate / fTrackEnergy;

    //rescale tracks
    for(Candidate *const &track : *fTowerTrackArray)
    {
      Candidate *new_track = static_cast<Candidate *>(track->Clone());
      new_track->AddCandidate(track);
      new_track->Momentum.SetPtEtaPhiM(new_track->Momentum.Pt() * rescaleFactor, new_track->Momentum.Eta(), new_track->Momentum.Phi(), new_track->Momentum.M());
      if(debug) cout << "  track Momentum: " << new_track->PID << ", " << new_track->Momentum.Pt() << ", " << new_track->Momentum.Eta() << ", " << new_track->Momentum.M() << endl;
      fEFlowTrackOutputArray->emplace_back(new_track);
    }
  }
}

//------------------------------------------------------------------------------

double DualReadoutCalorimeter::LogNormal(double mean, double sigma)
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

double DualReadoutCalorimeter::TruncatedGaussian(double mean, double sigma)
{
  double result = -1;
  if(mean > 0.0)
  {
    while(result < 0.0)
    {
      result = gRandom->Gaus(mean, sigma);
    }
    return result;
  }
  else
  {
    return 0.0;
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("DualReadoutCalorimeter", DualReadoutCalorimeter);
