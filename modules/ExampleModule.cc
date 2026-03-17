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
#include "classes/DelphesModuleFactory.h"

#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TRandom3.h>

#include <deque>

using namespace std;

class ExampleModule: public DelphesModule
{
public:
  ExampleModule() : fFormula(std::make_unique<DelphesFormula>()) {}

  void Init() override;
  void Process() override;
  void Finish() override {}

private:
  Int_t fIntParam;
  Double_t fDoubleParam;

  std::deque<Double_t> fArrayParam;

  const std::unique_ptr<DelphesFormula> fFormula; //!

  const TObjArray *fInputArray{nullptr}; //!
  std::unique_ptr<TIterator> fItInputArray; //!

  TObjArray *fOutputArray{nullptr}; //!
};

//------------------------------------------------------------------------------

void ExampleModule::Init()
{
  // read parameters

  fIntParam = GetInt("IntParam", 10);
  fDoubleParam = GetDouble("DoubleParam", 1.0);
  fFormula->Compile(GetString("EfficiencyFormula", "0.4"));

  ExRootConfParam param = GetParam("ArrayParam");
  Long_t i, size;
  fArrayParam.clear();

  size = param.GetSize();
  for(i = 0; i < size; ++i)
  {
    fArrayParam.push_back(param[i].GetDouble());
  }

  // import input array(s)

  fInputArray = ImportArray(GetString("InputArray", "FastJetFinder/jets"));
  fItInputArray.reset(fInputArray->MakeIterator());

  // create output array(s)

  fOutputArray = ExportArray(GetString("OutputArray", "jets"));
}

//------------------------------------------------------------------------------

void ExampleModule::Process()
{
  Candidate *candidate = nullptr;
  TLorentzVector candidatePosition, candidateMomentum;

  // loop over all input candidates
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    candidatePosition = candidate->Position;
    candidateMomentum = candidate->Momentum;

    // apply an efficency formula
    if(gRandom->Uniform() <= fFormula->Eval(candidateMomentum.Pt(), candidatePosition.Eta()))
    {
      fOutputArray->Add(candidate);
    }
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("ExampleModule", ExampleModule);
