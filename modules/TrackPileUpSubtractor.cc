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

#include <TLorentzVector.h>

using namespace std;

class TrackPileUpSubtractor: public DelphesModule
{
public:
  explicit TrackPileUpSubtractor(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fPTMin(Steer<double>("PTMin", 0.)),
    fFormula(std::make_unique<DelphesFormula>())
  {
    // read resolution formula in m
    fFormula->Compile(Steer<std::string>("ZVertexResolution", "0.001"));
  }

  void Init() override
  {
    fVertexInputArray = ImportArray(Steer<std::string>("VertexInputArray", "PileUpMerger/vertices"));
    // import arrays with output from other modules
    for(const std::pair<std::string, std::string> &inputArray :
      Steer<std::vector<std::pair<std::string, std::string> > >("InputArray"))
      fInputMap.emplace_back(std::make_pair(
        ImportArray(inputArray.first),
        ExportArray(inputArray.second)));
  }
  void Process() override;

private:
  const double fPTMin;
  const std::unique_ptr<DelphesFormula> fFormula; //!

  CandidatesCollection fVertexInputArray; //!
  std::vector<std::pair<CandidatesCollection, CandidatesCollection> > fInputMap; //!
};

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Process()
{
  for(const auto &[input_collection, output_collection] : fInputMap)
    output_collection->clear();

  double zvtx = 0.;
  // find z position of primary vertex
  for(Candidate *const &candidate : *fVertexInputArray)
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
    for(Candidate *const &candidate : *input_collection)
    {
      Candidate *particle = static_cast<Candidate *>(candidate->GetCandidates().at(0));
      const TLorentzVector &candidateMomentum = particle->Momentum;

      const double eta = candidateMomentum.Eta();
      const double pt = candidateMomentum.Pt();
      const double phi = candidateMomentum.Phi();
      const double e = candidateMomentum.E();

      const double z = particle->Position.Z();

      // apply pile-up subtraction
      // assume perfect pile-up subtraction for tracks outside fZVertexResolution

      if(candidate->Charge != 0 && candidate->IsPU && std::fabs(z - zvtx) > fFormula->Eval(pt, eta, phi, e) * 1.0e3)
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
