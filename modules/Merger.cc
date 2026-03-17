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

/** \class Merger
 *
 *  Merges multiple input arrays into one output array
 *  and sums transverse momenta of all input objects.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"
#include "classes/DelphesModuleFactory.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include <TLorentzVector.h>
#include <TObjArray.h>

#include <vector>

using namespace std;

class Merger: public DelphesModule
{
public:
  Merger() = default;

  void Init() override;
  void Process() override;
  void Finish() override;

private:
  std::vector<std::unique_ptr<TIterator> > fInputList; //!

  TObjArray *fOutputArray{nullptr}; //!
  TObjArray *fMomentumOutputArray{nullptr}; //!
  TObjArray *fEnergyOutputArray{nullptr}; //!
};

//------------------------------------------------------------------------------

void Merger::Init()
{
  // import arrays with output from other modules

  ExRootConfParam param = GetParam("InputArray");
  Long_t i, size;
  const TObjArray *array;

  size = param.GetSize();
  for(i = 0; i < size; ++i)
  {
    array = ImportArray(param[i].GetString());
    fInputList.push_back(std::unique_ptr<TIterator>(array->MakeIterator()));
  }

  // create output arrays

  fOutputArray = ExportArray(GetString("OutputArray", "candidates"));

  fMomentumOutputArray = ExportArray(GetString("MomentumOutputArray", "momentum"));

  fEnergyOutputArray = ExportArray(GetString("EnergyOutputArray", "energy"));
}

//------------------------------------------------------------------------------

void Merger::Finish()
{
  fInputList.clear();
}

//------------------------------------------------------------------------------

void Merger::Process()
{
  Candidate *candidate = nullptr;
  TLorentzVector momentum;
  Double_t sumPT, sumE;
  vector<std::unique_ptr<TIterator> >::iterator itInputList;

  DelphesFactory *factory = GetFactory();

  momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  sumPT = 0;
  sumE = 0;

  // loop over all input arrays
  for(itInputList = fInputList.begin(); itInputList != fInputList.end(); ++itInputList)
  {
    auto &iterator = *itInputList;

    // loop over all candidates
    iterator->Reset();
    while((candidate = static_cast<Candidate *>(iterator->Next())))
    {
      const TLorentzVector &candidateMomentum = candidate->Momentum;

      momentum += candidateMomentum;
      sumPT += candidateMomentum.Pt();
      sumE += candidateMomentum.E();

      fOutputArray->Add(candidate);
    }
  }

  candidate = factory->NewCandidate();

  candidate->Position.SetXYZT(0.0, 0.0, 0.0, 0.0);
  candidate->Momentum = momentum;

  fMomentumOutputArray->Add(candidate);

  candidate = factory->NewCandidate();

  candidate->Position.SetXYZT(0.0, 0.0, 0.0, 0.0);
  candidate->Momentum.SetPtEtaPhiE(sumPT, 0.0, 0.0, sumE);

  fEnergyOutputArray->Add(candidate);
}

//------------------------------------------------------------------------------

REGISTER_MODULE("Merger", Merger);
