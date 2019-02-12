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

#ifndef DenseTrackFilter_h
#define DenseTrackFilter_h

/** \class DenseTrackFilter
 *
 *  Applies reconstruction inefficiency on tracks in dense environment
 *
 *  \author M. Selvaggi - CERN
 *
 */

#include "classes/DelphesModule.h"

#include <map>
#include <set>
#include <vector>

class TObjArray;
class DelphesFormula;
class Candidate;

class DenseTrackFilter: public DelphesModule
{
public:
  DenseTrackFilter();
  ~DenseTrackFilter();

  void Init();
  void Process();
  void Finish();

private:
  typedef std::map<Double_t, std::set<Double_t> > TBinMap; //!

  Candidate *fBestTrack;

  Int_t fTowerTrackHits;

  Double_t fEtaPhiRes;

  TBinMap fBinMap; //!

  std::vector<Double_t> fEtaBins;
  std::vector<std::vector<Double_t> *> fPhiBins;

  std::vector<Long64_t> fTowerHits;

  TIterator *fItTrackInputArray; //!

  const TObjArray *fTrackInputArray; //!
  TObjArray *fTrackOutputArray; //!

  TObjArray *fChargedHadronOutputArray; //!
  TObjArray *fElectronOutputArray; //!
  TObjArray *fMuonOutputArray; //!

  void FillTrack();
  ClassDef(DenseTrackFilter, 1)
};

#endif
