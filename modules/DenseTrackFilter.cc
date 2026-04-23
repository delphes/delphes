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

/** \class DenseTrackFilter
 *
 *  Applies reconstruction inefficiency on tracks in dense environment
 *
 *  \author M. Selvaggi - CERN
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

#include <set>

using namespace std;

class DenseTrackFilter: public DelphesModule
{
public:
  explicit DenseTrackFilter(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fEtaPhiRes(Steer<double>("EtaPhiRes", 0.003)), // Eta x Phi smearing to be applied
    fPhiBins(Steer<std::unordered_map<double, std::vector<double> > >("EtaPhiBins"))
  {
    for(const std::pair<double, std::vector<double> > etaPhiBins : fPhiBins)
      fEtaBins.insert(etaPhiBins.first);
  }

  void Init() override
  {
    fTrackInputArray = ImportArray(Steer<std::string>("TrackInputArray", "TrackMergerProp/tracks"));
    fTrackOutputArray = ExportArray(Steer<std::string>("TrackOutputArray", "tracks"));
    fChargedHadronOutputArray = ExportArray(Steer<std::string>("ChargedHadronOutputArray", "chargedHadrons"));
    fElectronOutputArray = ExportArray(Steer<std::string>("ElectronOutputArray", "electrons"));
    fMuonOutputArray = ExportArray(Steer<std::string>("MuonOutputArray", "muons"));
  }
  void Process() override;

private:
  void FillTrack();

  const double fEtaPhiRes;

  typedef std::map<double, std::set<double> > TBinMap; //!

  Candidate *fBestTrack{nullptr};

  int fTowerTrackHits;

  TBinMap fBinMap; //!

  const std::unordered_map<double, std::vector<double> > fPhiBins;
  std::set<double> fEtaBins;

  std::vector<long long> fTowerHits;

  CandidatesCollection fTrackInputArray; //!

  CandidatesCollection fTrackOutputArray; //!
  CandidatesCollection fChargedHadronOutputArray; //!
  CandidatesCollection fElectronOutputArray; //!
  CandidatesCollection fMuonOutputArray; //!
};

//------------------------------------------------------------------------------

void DenseTrackFilter::Process()
{
  fTrackOutputArray->clear();
  fChargedHadronOutputArray->clear();
  fElectronOutputArray->clear();
  fMuonOutputArray->clear();

  fTowerHits.clear();

  // loop over all tracks
  size_t number = 0;
  for(Candidate *const &track : *fTrackInputArray)
  {
    const TLorentzVector &trackPosition = track->Position;

    // find eta bin [1, fEtaBins.size - 1]
    std::set<double>::iterator itEtaBin = std::lower_bound(fEtaBins.begin(), fEtaBins.end(), trackPosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    const short etaBin = std::distance(fEtaBins.begin(), itEtaBin);

    // phi bins for given eta bin
    const std::vector<double> &phiBins = fPhiBins.at(*itEtaBin);

    // find phi bin [1, phiBins.size - 1]
    std::vector<double>::const_iterator itPhiBin = std::lower_bound(phiBins.begin(), phiBins.end(), trackPosition.Phi());
    if(itPhiBin == phiBins.begin() || itPhiBin == phiBins.end()) continue; // outside range
    const short phiBin = std::distance(phiBins.begin(), itPhiBin);

    const short flags = 1;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for track number}
    const unsigned long long towerHit = (Long64_t(etaBin) << 48) | (Long64_t(phiBin) << 32) | (Long64_t(flags) << 24) | Long64_t(number);

    fTowerHits.push_back(towerHit);
    ++number;
  }

  // all hits are sorted first by eta bin number, then by phi bin number,
  // then by flags and then by particle or track number
  std::sort(fTowerHits.begin(), fTowerHits.end());

  // loop over all hits
  unsigned long long towerEtaPhi = 0;
  fBestTrack = 0;
  double ptmax = 0.;
  fTowerTrackHits = 0;

  for(const long long &towerHit : fTowerHits)
  {
    const short flags = (towerHit >> 24) & 0x00000000000000FFLL;
    const size_t number = (towerHit) & 0x0000000000FFFFFFLL;
    const unsigned long long hitEtaPhi = towerHit >> 32;

    if(towerEtaPhi != hitEtaPhi)
    {
      // switch to next tower
      towerEtaPhi = hitEtaPhi;

      // saving track with highest pT that hit the tower
      FillTrack();

      ptmax = 0.;
      fTowerTrackHits = 0;
      fBestTrack = 0;
    }
    // check for track hits

    if(flags & 1)
    {
      ++fTowerTrackHits;
      Candidate *track = static_cast<Candidate *>(fTrackInputArray->at(number));
      const TLorentzVector &momentum = track->Momentum;

      if(momentum.Pt() > ptmax)
      {
        ptmax = momentum.Pt();
        fBestTrack = track;
      }
      continue;
    }
  }

  // here fill last tower
  FillTrack();
}

//------------------------------------------------------------------------------

void DenseTrackFilter::FillTrack()
{
  // saving track with highest pT that hit the tower
  if(fTowerTrackHits < 1) return;

  const size_t numberOfCandidates = fBestTrack->GetCandidates().size();
  if(numberOfCandidates < 1) return;

  Candidate *track = static_cast<Candidate *>(fBestTrack->GetCandidates().at(numberOfCandidates - 1));
  Candidate *new_candidate = static_cast<Candidate *>(track->Clone());

  const double pt = new_candidate->Momentum.Pt();
  const double eta = gRandom->Gaus(new_candidate->Momentum.Eta(), fEtaPhiRes);
  const double phi = gRandom->Gaus(new_candidate->Momentum.Phi(), fEtaPhiRes);
  const double m = new_candidate->Momentum.M();
  new_candidate->Momentum.SetPtEtaPhiM(pt, eta, phi, m);
  new_candidate->AddCandidate(track);

  fTrackOutputArray->emplace_back(new_candidate);
  switch(std::abs(new_candidate->PID))
  {
  case 11:
    fElectronOutputArray->emplace_back(new_candidate);
    break;
  case 13:
    fMuonOutputArray->emplace_back(new_candidate);
    break;
  default:
    fChargedHadronOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("DenseTrackFilter", DenseTrackFilter);
