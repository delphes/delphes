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

/** \class EnergyScale
 *
 *  Applies energy scale.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>

using namespace std;

class EnergyScale: public DelphesModule
{
public:
  explicit EnergyScale(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fFormula(std::make_unique<DelphesFormula>())
  {
    // read resolution formula
    fFormula->Compile(Steer<std::string>("ScaleFormula", "0.0"));
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "FastJetFinder/jets")); // import input arrays
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "jets")); // create output array
  }
  void Process() override
  {
    fOutputArray->clear();
    for(Candidate *const &candidate : *fInputArray)
    {
      TLorentzVector momentum = candidate->Momentum;

      const double scale = fFormula->Eval(momentum.Pt(), momentum.Eta(), momentum.Phi(), momentum.E());
      if(scale > 0.0) momentum *= scale;

      Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());
      new_candidate->Momentum = momentum;

      fOutputArray->emplace_back(new_candidate);
    }
  }

private:
  const std::unique_ptr<DelphesFormula> fFormula; //!

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

REGISTER_MODULE("EnergyScale", EnergyScale);
