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

/** \class LeptonDressing
 *
 *  \author P. Demin && A. Mertens - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>

using namespace std;

class LeptonDressing: public DelphesModule
{
public:
  LeptonDressing() = default;

  void Init() override;
  void Process() override;

private:
  Double_t fDeltaR;

  CandidatesCollection fDressingInputArray; //!
  CandidatesCollection fCandidateInputArray; //!

  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void LeptonDressing::Init()
{
  fDeltaR = GetDouble("DeltaRMax", 0.4);

  // import input arrays
  fDressingInputArray = ImportArray(GetString("DressingInputArray", "Calorimeter/photons"));
  fCandidateInputArray = ImportArray(GetString("CandidateInputArray", "UniqueObjectFinder/electrons"));

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "electrons"));
}

//------------------------------------------------------------------------------

void LeptonDressing::Process()
{
  fOutputArray->clear();

  TLorentzVector momentum;

  // loop over all input candidate
  for(const auto &candidate : *fCandidateInputArray)
  {
    const TLorentzVector &candidateMomentum = candidate->Momentum;

    // loop over all input tracks
    momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    for(const auto &dressing : *fDressingInputArray)
    {
      const TLorentzVector &dressingMomentum = dressing->Momentum;
      if(dressingMomentum.Pt() > 0.1)
      {
        if(candidateMomentum.DeltaR(dressingMomentum) <= fDeltaR)
        {
          momentum += dressingMomentum;
        }
      }
    }

    auto *new_candidate = static_cast<Candidate *>(candidate->Clone());

    new_candidate->Momentum += momentum;
    new_candidate->AddCandidate(candidate);

    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("LeptonDressing", LeptonDressing);
