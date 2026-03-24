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

/** \class CscClusterEfficiency
 *
 *  This module is specific to the CMS paper searching for neutral LLPs in the CMS endcap muon detectors: https://arxiv.org/abs/2107.04838
 *  It is implemented based on the ClusterEfficiency parameterization function provided in the HEPData entry of the paper: https://www.hepdata.net/record/104408
 *
 *  \author Christina Wang
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesCscClusterFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

using namespace std;

class CscClusterEfficiency: public DelphesModule
{
public:
  explicit CscClusterEfficiency(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fFormula(std::make_unique<DelphesCscClusterFormula>())
  {
    // read CscClusterEfficiency formula
    fFormula->Compile(Steer<std::string>("EfficiencyFormula", "1.0"));
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "ParticlePropagator/stableParticles")); // import input array
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "stableParticles")); // create output array
  }
  void Process() override;

private:
  const std::unique_ptr<DelphesCscClusterFormula> fFormula; //!

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void CscClusterEfficiency::Process()
{
  fOutputArray->clear();
  for(Candidate *const &candidate : *fInputArray)
  {
    const TLorentzVector &candidateDecayPosition = candidate->DecayPosition;
    double decayZ = abs(candidateDecayPosition.Z());
    double decayR = sqrt(pow(candidateDecayPosition.X(), 2) + pow(candidateDecayPosition.Y(), 2));
    double Ehad = candidate->Ehad;
    double Eem = candidate->Eem;
    // apply an efficency formula
    if(gRandom->Uniform() > fFormula->Eval(decayR, decayZ, Ehad, Eem)) continue;

    fOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("CscClusterEfficiency", CscClusterEfficiency);
