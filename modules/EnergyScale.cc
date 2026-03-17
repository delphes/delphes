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
#include "classes/DelphesModuleFactory.h"

#include <TLorentzVector.h>
#include <TObjArray.h>

using namespace std;

class EnergyScale: public DelphesModule
{
public:
  EnergyScale() : fFormula(std::make_unique<DelphesFormula>()) {}

  void Init() override;
  void Process() override;

private:
  const std::unique_ptr<DelphesFormula> fFormula; //!

  const TObjArray *fInputArray{nullptr}; //!
  std::unique_ptr<TIterator> fItInputArray; //!

  TObjArray *fOutputArray{nullptr}; //!
};

//------------------------------------------------------------------------------

void EnergyScale::Init()
{
  // read resolution formula
  fFormula->Compile(GetString("ScaleFormula", "0.0"));

  // import input arrays
  fInputArray = ImportArray(GetString("InputArray", "FastJetFinder/jets"));
  fItInputArray.reset(fInputArray->MakeIterator());

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "jets"));
}

//------------------------------------------------------------------------------

void EnergyScale::Process()
{
  Candidate *candidate = nullptr;
  TLorentzVector momentum;
  Double_t scale;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    momentum = candidate->Momentum;

    scale = fFormula->Eval(momentum.Pt(), momentum.Eta(), momentum.Phi(), momentum.E());

    if(scale > 0.0) momentum *= scale;

    candidate = static_cast<Candidate *>(candidate->Clone());
    candidate->Momentum = momentum;

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("EnergyScale", EnergyScale);
