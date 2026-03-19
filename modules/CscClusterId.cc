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

/** \class CscClusterId
 *
 *  This module is specific to the CMS paper searching for neutral LLPs in the CMS endcap muon detectors: https://arxiv.org/abs/2107.04838
 *  It is implemented based on the cut_based_id.py function provided in the HEPData entry of the paper: https://www.hepdata.net/record/104408
 *
 *  \author Christina Wang
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesCscClusterFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>

using namespace std;

class CscClusterId: public DelphesModule
{
public:
  CscClusterId() :
    fFormula(std::make_unique<DelphesCscClusterFormula>()),
    fEtaFormula(std::make_unique<DelphesCscClusterFormula>()) {}

  void Init() override;
  void Process() override;

private:
  const std::unique_ptr<DelphesCscClusterFormula> fFormula; //!
  const std::unique_ptr<DelphesCscClusterFormula> fEtaFormula; //!
  Double_t fEtaCutMax;

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void CscClusterId::Init()
{
  // read efficiency formula
  fFormula->Compile(GetString("EfficiencyFormula", "1.0"));
  fEtaFormula->Compile(GetString("EtaCutFormula", "1.0"));
  fEtaCutMax = GetDouble("EtaCutMax", 999.0);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void CscClusterId::Process()
{
  fOutputArray->clear();

  Double_t Ehad, decayR, decayZ, NStationEff, eta;
  Double_t signPz, cosTheta;

  for(Candidate *const &candidate : *fInputArray)
  {
    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &candidateDecayPosition = candidate->DecayPosition;
    decayZ = abs(candidateDecayPosition.Z());
    decayR = sqrt(pow(candidateDecayPosition.X(), 2) + pow(candidateDecayPosition.Y(), 2));
    Ehad = candidate->Ehad;

    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signPz * 999.9 : momentum.Eta());

    // calculate the NStation > 1 efficiency, implemented according to Additional Figure 8 in HEPData
    NStationEff = fFormula->Eval(decayR, decayZ, Ehad);

    // depending on the decay region (station Number), different eta cut is applied, implemented based on cut_based_id.py in HEPData
    float eta_cut = fEtaFormula->Eval(decayR, decayZ);
    if(gRandom->Uniform() > NStationEff * (abs(eta) < fEtaCutMax) + (1.0 - NStationEff) * (abs(eta) < eta_cut)) continue;

    fOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("CscClusterId", CscClusterId);
