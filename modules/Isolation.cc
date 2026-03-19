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
#include <TMath.h>

using namespace std;

//------------------------------------------------------------------------------

class IsolationClassifier: public ExRootClassifier
{
public:
  IsolationClassifier() {}

  Int_t GetCategory(TObject *object)
  {
    Candidate *track = static_cast<Candidate *>(object);
    if(const TLorentzVector &momentum = track->Momentum; momentum.Pt() < fPTMin) return -1;
    return 0;
  }

  Double_t fPTMin;
};

//------------------------------------------------------------------------------

class Isolation: public DelphesModule
{
public:
  Isolation() : fClassifier(std::make_unique<IsolationClassifier>()) {}

  void Init() override;
  void Process() override;

private:
  Double_t fDeltaRMax;
  Double_t fPTRatioMax;
  Double_t fPTSumMax;
  Double_t fDeltaRMin;
  Bool_t fUsePTSum;
  Bool_t fUseRhoCorrection;
  Bool_t fUseMiniCone;

  const std::unique_ptr<IsolationClassifier> fClassifier; //!

  std::unique_ptr<DelphesFilter> fFilter;

  CandidatesCollection fIsolationInputArray; //!
  CandidatesCollection fCandidateInputArray; //!
  CandidatesCollection fRhoInputArray; //!

  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void Isolation::Init()
{
  const char *rhoInputArrayName;

  fDeltaRMax = GetDouble("DeltaRMax", 0.5);
  fPTRatioMax = GetDouble("PTRatioMax", 0.1);
  fPTSumMax = GetDouble("PTSumMax", 5.0);
  fUsePTSum = GetBool("UsePTSum", false);
  fUseRhoCorrection = GetBool("UseRhoCorrection", true);
  fDeltaRMin = GetDouble("DeltaRMin", 0.01);
  fUseMiniCone = GetBool("UseMiniCone", false);
  fClassifier->fPTMin = GetDouble("PTMin", 0.5);

  // import input arrays
  fIsolationInputArray = ImportArray(GetString("IsolationInputArray", "Delphes/partons"));
  fFilter = std::make_unique<DelphesFilter>(fIsolationInputArray);

  fCandidateInputArray = ImportArray(GetString("CandidateInputArray", "Calorimeter/electrons"));

  rhoInputArrayName = GetString("RhoInputArray", "");
  if(rhoInputArrayName[0] != '\0')
    fRhoInputArray = ImportArray(rhoInputArrayName);

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "electrons"));
}

//------------------------------------------------------------------------------

void Isolation::Process()
{
  fOutputArray->clear();

  CandidatesCollection isolationArray;
  Double_t sumChargedNoPU, sumChargedPU, sumNeutral, sumAllParticles;
  Double_t sumDBeta, ratioDBeta, sumRhoCorr, ratioRhoCorr, sum, ratio;
  Bool_t pass = kFALSE;
  Double_t eta = 0.0;
  Double_t rho = 0.0;

  // select isolation objects
  fFilter->Reset();
  isolationArray = fFilter->GetSubArray(fClassifier.get(), 0);

  // loop over all input jets
  for(Candidate *const &candidate : *fCandidateInputArray)
  {
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = TMath::Abs(candidateMomentum.Eta());

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
    sumDBeta = sumChargedNoPU + TMath::Max(sumNeutral - 0.5 * sumChargedPU, 0.0);
    sumRhoCorr = sumChargedNoPU + TMath::Max(sumNeutral - TMath::Max(rho, 0.0) * fDeltaRMax * fDeltaRMax * TMath::Pi(), 0.0);
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
