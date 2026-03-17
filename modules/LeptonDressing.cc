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
#include "classes/DelphesModuleFactory.h"

#include <TLorentzVector.h>
#include <TObjArray.h>

using namespace std;

class LeptonDressing: public DelphesModule
{
public:
  LeptonDressing() = default;

  void Init() override;
  void Process() override;

private:
  Double_t fDeltaR;

  const TObjArray *fDressingInputArray{nullptr}; //!
  std::unique_ptr<TIterator> fItDressingInputArray; //!

  const TObjArray *fCandidateInputArray{nullptr}; //!
  std::unique_ptr<TIterator> fItCandidateInputArray; //!

  TObjArray *fOutputArray{nullptr}; //!
};

//------------------------------------------------------------------------------

void LeptonDressing::Init()
{
  fDeltaR = GetDouble("DeltaRMax", 0.4);

  // import input arrasy
  fDressingInputArray = ImportArray(GetString("DressingInputArray", "Calorimeter/photons"));
  fItDressingInputArray.reset(fDressingInputArray->MakeIterator());

  fCandidateInputArray = ImportArray(GetString("CandidateInputArray", "UniqueObjectFinder/electrons"));
  fItCandidateInputArray.reset(fCandidateInputArray->MakeIterator());

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "electrons"));
}

//------------------------------------------------------------------------------

void LeptonDressing::Process()
{
  Candidate *candidate = nullptr, *dressing = nullptr, *mother = nullptr;
  TLorentzVector momentum;

  // loop over all input candidate
  fItCandidateInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItCandidateInputArray->Next())))
  {
    const TLorentzVector &candidateMomentum = candidate->Momentum;

    // loop over all input tracks
    fItDressingInputArray->Reset();
    momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    while((dressing = static_cast<Candidate *>(fItDressingInputArray->Next())))
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

    mother = candidate;
    candidate = static_cast<Candidate *>(candidate->Clone());

    candidate->Momentum += momentum;
    candidate->AddCandidate(mother);

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("LeptonDressing", LeptonDressing);
