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

/** \class TimeSmearing
 *
 *  Performs time smearing.
 *
 *  \author Michele Selvaggi - CERN
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

using namespace std;

class TimeSmearing: public DelphesModule
{
public:
  explicit TimeSmearing(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fResolutionFormula(std::make_unique<DelphesFormula>())
  {
    // read time resolution formula in seconds
    fResolutionFormula->Compile(Steer<std::string>("TimeResolution", "30e-12"));
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "MuonMomentumSmearing/muons"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "tracks"));
  }
  void Process() override;

private:
  const std::unique_ptr<DelphesFormula> fResolutionFormula;

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void TimeSmearing::Process()
{
  fOutputArray->clear();

  const double c_light = 2.99792458E8;

  for(Candidate *const &candidate : *fInputArray)
  {
    const TLorentzVector &candidateFinalPosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;

    const double tf = candidateFinalPosition.T() * 1.0E-3 / c_light;

    const double eta = candidateMomentum.Eta();
    const double energy = candidateMomentum.E();

    // apply smearing formula
    const double timeResolution = fResolutionFormula->Eval(0.0, eta, 0.0, energy);
    const double tf_smeared = gRandom->Gaus(tf, timeResolution);

    Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());

    new_candidate->Position.SetT(tf_smeared * 1.0E3 * c_light);
    new_candidate->ErrorT = timeResolution * 1.0E3 * c_light;

    new_candidate->AddCandidate(candidate);
    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("TimeSmearing", TimeSmearing);
