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

/** \class ExampleModule
 *
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

#include <deque>

using namespace std;

class ExampleModule: public DelphesModule
{
public:
  explicit ExampleModule(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fIntParam(Steer<int>("IntParam", 10)),
    fDoubleParam(Steer<double>("DoubleParam", 1.0)),
    fFormula(std::make_unique<DelphesFormula>())
  {
    fFormula->Compile(Steer<std::string>("EfficiencyFormula", "0.4"));
    for(const double &param : Steer<std::vector<double> >("ArrayParam"))
      fArrayParam.push_back(param);
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "FastJetFinder/jets")); // import input array(s)
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "jets")); // create output array(s)
  }
  void Process() override;
  void Finish() override {}

private:
  const int fIntParam;
  const double fDoubleParam;

  const std::unique_ptr<DelphesFormula> fFormula; //!

  std::deque<double> fArrayParam;

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void ExampleModule::Process()
{
  fOutputArray->clear();

  // loop over all input candidates
  for(Candidate *const &candidate : *fInputArray)
  {
    // apply an efficency formula
    if(gRandom->Uniform() <= fFormula->Eval(candidate->Momentum.Pt(), candidate->Position.Eta()))
      fOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("ExampleModule", ExampleModule);
