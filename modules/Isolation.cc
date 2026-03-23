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

/** \class Isolation
 *
 *  Sums transverse momenta of isolation objects (tracks, calorimeter towers, etc)
 *  within a DeltaR cone around a candidate and calculates fraction of this sum
 *  to the candidate's transverse momentum. outputs candidates that have
 *  the transverse momenta fraction within (PTRatioMin, PTRatioMax].
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFilter.h"
#include "classes/DelphesModule.h"

#include "ExRootAnalysis/ExRootClassifier.h"

#include <TLorentzVector.h>

using namespace std;

//------------------------------------------------------------------------------

class IsolationClassifier: public ExRootClassifier
{
public:
  IsolationClassifier() {}

  int GetCategory(TObject *object)
  {
    Candidate *track = static_cast<Candidate *>(object);
    if(const TLorentzVector &momentum = track->Momentum; momentum.Pt() < fPTMin) return -1;
    return 0;
  }

  double fPTMin;
};

//------------------------------------------------------------------------------

class Isolation: public DelphesModule
{
public:
  explicit Isolation(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fDeltaRMax(Steer<double>("DeltaRMax", 0.5)),
    fPTRatioMax(Steer<double>("PTRatioMax", 0.1)),
    fPTSumMax(Steer<double>("PTSumMax", 5.0)),
    fDeltaRMin(Steer<double>("DeltaRMin", 0.01)),
    fUsePTSum(Steer<bool>("UsePTSum", false)),
    fUseRhoCorrection(Steer<bool>("UseRhoCorrection", true)),
    fUseMiniCone(Steer<bool>("UseMiniCone", false)),
    fClassifier(std::make_unique<IsolationClassifier>())
  {
    fClassifier->fPTMin = Steer<double>("PTMin", 0.5);
  }

  void Init() override
  {
    fIsolationInputArray = ImportArray(Steer<std::string>("IsolationInputArray", "Delphes/partons"));
    fCandidateInputArray = ImportArray(Steer<std::string>("CandidateInputArray", "Calorimeter/electrons"));
    if(const std::string rhoInputArrayName = Steer<std::string>("RhoInputArray", ""); !rhoInputArrayName.empty())
      fRhoInputArray = ImportArray(rhoInputArrayName);
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "electrons"));

    fFilter = std::make_unique<DelphesFilter>(fIsolationInputArray);
  }
  void Process() override;

private:
  const double fDeltaRMax;
  const double fPTRatioMax;
  const double fPTSumMax;
  const double fDeltaRMin;
  const bool fUsePTSum;
  const bool fUseRhoCorrection;
  const bool fUseMiniCone;

  const std::unique_ptr<IsolationClassifier> fClassifier; //!

  std::unique_ptr<DelphesFilter> fFilter;

  CandidatesCollection fIsolationInputArray; //!
  CandidatesCollection fCandidateInputArray; //!
  CandidatesCollection fRhoInputArray; //!

  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void Isolation::Process()
{
  fOutputArray->clear();

  CandidatesCollection isolationArray;
  double sumChargedNoPU, sumChargedPU, sumNeutral, sumAllParticles;
  double sumDBeta, ratioDBeta, sumRhoCorr, ratioRhoCorr, sum, ratio;
  bool pass = kFALSE;
  double eta = 0.0;
  double rho = 0.0;

  // select isolation objects
  fFilter->Reset();
  isolationArray = fFilter->GetSubArray(fClassifier.get(), 0);

  // loop over all input jets
  for(Candidate *const &candidate : *fCandidateInputArray)
  {
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = std::fabs(candidateMomentum.Eta());

    // find rho
    rho = 0.0;
    if(fRhoInputArray)
    {
      for(Candidate *const &object : *fRhoInputArray)
        if(eta >= object->Edges[0] && eta < object->Edges[1])
          rho = object->Momentum.Pt();
    }

    // loop over all input tracks

    sumNeutral = 0.0;
    sumChargedNoPU = 0.0;
    sumChargedPU = 0.0;
    sumAllParticles = 0.0;

    for(Candidate *const &isolation : *isolationArray)
    {
      const TLorentzVector &isolationMomentum = isolation->Momentum;

      if(fUseMiniCone)
        pass = candidateMomentum.DeltaR(isolationMomentum) <= fDeltaRMax && candidateMomentum.DeltaR(isolationMomentum) > fDeltaRMin;
      else
        pass = candidateMomentum.DeltaR(isolationMomentum) <= fDeltaRMax && candidate->GetUniqueID() != isolation->GetUniqueID();

      if(pass)
      {
        sumAllParticles += isolationMomentum.Pt();
        if(isolation->Charge != 0)
        {
          if(isolation->IsRecoPU)
            sumChargedPU += isolationMomentum.Pt();
          else
            sumChargedNoPU += isolationMomentum.Pt();
        }
        else
          sumNeutral += isolationMomentum.Pt();
      }
    }

    // find rho
    rho = 0.0;
    if(fRhoInputArray)
    {
      for(Candidate *const &object : *fRhoInputArray)
        if(eta >= object->Edges[0] && eta < object->Edges[1])
          rho = object->Momentum.Pt();
    }

    // correct sum for pile-up contamination
    sumDBeta = sumChargedNoPU + std::max(sumNeutral - 0.5 * sumChargedPU, 0.0);
    sumRhoCorr = sumChargedNoPU + std::max(sumNeutral - std::max(rho, 0.0) * fDeltaRMax * fDeltaRMax * M_PI, 0.0);
    ratioDBeta = sumDBeta / candidateMomentum.Pt();
    ratioRhoCorr = sumRhoCorr / candidateMomentum.Pt();

    candidate->IsolationVar = ratioDBeta;
    candidate->IsolationVarRhoCorr = ratioRhoCorr;
    candidate->SumPtCharged = sumChargedNoPU;
    candidate->SumPtNeutral = sumNeutral;
    candidate->SumPtChargedPU = sumChargedPU;
    candidate->SumPt = sumAllParticles;

    sum = fUseRhoCorrection ? sumRhoCorr : sumDBeta;
    if(fUsePTSum && sum > fPTSumMax) continue;

    ratio = fUseRhoCorrection ? ratioRhoCorr : ratioDBeta;
    if(!fUsePTSum && ratio > fPTRatioMax) continue;

    fOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("Isolation", Isolation);
