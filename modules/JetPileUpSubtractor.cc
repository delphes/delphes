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

/** \class JetPileUpSubtractor
 *
 *  Subtract pile-up contribution from jets using the fastjet area method
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"
#include "classes/DelphesModuleFactory.h"

#include <TLorentzVector.h>
#include <TObjArray.h>

using namespace std;

class JetPileUpSubtractor: public DelphesModule
{
public:
  JetPileUpSubtractor() = default;

  void Init() override;
  void Process() override;

private:
  Double_t fJetPTMin;

  const TObjArray *fJetInputArray{nullptr}; //!
  std::unique_ptr<TIterator> fItJetInputArray; //!

  const TObjArray *fRhoInputArray{nullptr}; //!
  std::unique_ptr<TIterator> fItRhoInputArray; //!

  TObjArray *fOutputArray{nullptr}; //!
};

//------------------------------------------------------------------------------

void JetPileUpSubtractor::Init()
{
  fJetPTMin = GetDouble("JetPTMin", 20.0);

  // import input arrays
  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray.reset(fJetInputArray->MakeIterator());

  fRhoInputArray = ImportArray(GetString("RhoInputArray", "Rho/rho"));
  fItRhoInputArray.reset(fRhoInputArray->MakeIterator());

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "jets"));
}

//------------------------------------------------------------------------------

void JetPileUpSubtractor::Process()
{
  Candidate *candidate = nullptr, *object = nullptr;
  TLorentzVector momentum, area;
  Double_t eta = 0.0;
  Double_t rho = 0.0;

  // loop over all input candidates
  fItJetInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItJetInputArray->Next())))
  {
    momentum = candidate->Momentum;
    area = candidate->Area;
    eta = momentum.Eta();

    // find rho
    rho = 0.0;
    if(fRhoInputArray)
    {
      fItRhoInputArray->Reset();
      while((object = static_cast<Candidate *>(fItRhoInputArray->Next())))
      {
        if(eta >= object->Edges[0] && eta < object->Edges[1])
        {
          rho = object->Momentum.Pt();
        }
      }
    }

    // apply pile-up correction
    if(momentum.Pt() <= rho * area.Pt()) continue;

    momentum -= rho * area;

    if(momentum.Pt() <= fJetPTMin) continue;

    candidate = static_cast<Candidate *>(candidate->Clone());
    candidate->Momentum = momentum;

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("JetPileUpSubtractor", JetPileUpSubtractor);
