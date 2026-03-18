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

/** \class TrackPileUpSubtractor
 *
 *  Subtract pile-up contribution from tracks.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"
#include "classes/DelphesModuleFactory.h"

#include <TLorentzVector.h>
#include <TMath.h>

using namespace std;

class TrackPileUpSubtractor: public DelphesModule
{
public:
  TrackPileUpSubtractor() : fFormula(std::make_unique<DelphesFormula>()) {}

  void Init() override;
  void Process() override;
  void Finish() override;

private:
  const std::unique_ptr<DelphesFormula> fFormula; //!

  Double_t fPTMin;

  std::vector<std::pair<CandidatesCollection, CandidatesCollection> > fInputMap; //!

  CandidatesCollection fVertexInputArray{nullptr}; //!
};

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Init()
{
  // import input array
  fVertexInputArray = ImportArray(GetString("VertexInputArray", "PileUpMerger/vertices"));

  // read resolution formula in m
  fFormula->Compile(GetString("ZVertexResolution", "0.001"));

  fPTMin = GetDouble("PTMin", 0.);

  // import arrays with output from other modules
  ExRootConfParam param = GetParam("InputArray");
  Long_t i, size;

  size = param.GetSize();
  for(i = 0; i < size / 2; ++i)
    fInputMap.emplace_back(std::make_pair(ImportArray(param[i * 2].GetString()), ExportArray(param[i * 2 + 1].GetString())));
}

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Finish()
{
  fInputMap.clear();
}

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Process()
{
  for(const auto &[input_collection, output_collection] : fInputMap)
    output_collection->clear();

  Double_t z, zvtx = 0;
  Double_t pt, eta, phi, e;

  // find z position of primary vertex
  for(const auto &candidate : *fVertexInputArray)
  {
    if(!candidate->IsPU)
    {
      zvtx = candidate->Position.Z();
      // break;
    }
  }

  // loop over all input arrays
  for(const auto &[input_collection, output_collection] : fInputMap)
  {
    // loop over all candidates
    for(const auto &candidate : *input_collection)
    {
      auto *particle = static_cast<Candidate *>(candidate->GetCandidates().at(0));
      const TLorentzVector &candidateMomentum = particle->Momentum;

      eta = candidateMomentum.Eta();
      pt = candidateMomentum.Pt();
      phi = candidateMomentum.Phi();
      e = candidateMomentum.E();

      z = particle->Position.Z();

      // apply pile-up subtraction
      // assume perfect pile-up subtraction for tracks outside fZVertexResolution

      if(candidate->Charge != 0 && candidate->IsPU && TMath::Abs(z - zvtx) > fFormula->Eval(pt, eta, phi, e) * 1.0e3)
      {
        candidate->IsRecoPU = 1;
      }
      else
      {
        candidate->IsRecoPU = 0;
        if(candidate->Momentum.Pt() > fPTMin) output_collection->emplace_back(candidate);
      }
    }
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("TrackPileUpSubtractor", TrackPileUpSubtractor);
