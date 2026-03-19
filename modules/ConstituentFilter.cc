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

/** \class ConstituentFilter
 *
 *  Drops all input objects that are not constituents of any jet.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

using namespace std;

class ConstituentFilter: public DelphesModule
{
public:
  ConstituentFilter() = default;

  void Init() override;
  void Process() override;

private:
  Double_t fJetPTMin;

  std::vector<CandidatesCollection> fInputList; //!
  std::vector<std::pair<CandidatesCollection, CandidatesCollection> > fInputMap; //!
};

//------------------------------------------------------------------------------

void ConstituentFilter::Init()
{
  ExRootConfParam param;
  Long_t i, size;

  fJetPTMin = GetDouble("JetPTMin", 0.0);

  // import input arrays
  param = GetParam("JetInputArray");
  size = param.GetSize();
  for(i = 0; i < size; ++i)
    fInputList.emplace_back(ImportArray(param[i].GetString()));

  param = GetParam("ConstituentInputArray");
  size = param.GetSize();
  for(i = 0; i < size / 2; ++i)
    fInputMap.emplace_back(std::make_pair(ImportArray(param[i * 2].GetString()), ExportArray(param[i * 2 + 1].GetString())));
}

//------------------------------------------------------------------------------

void ConstituentFilter::Process()
{
  for(auto &[input_collection, output_collection] : fInputMap) output_collection->clear();

  // loop over all jet input arrays
  for(const CandidatesCollection &input_collection : fInputList)
  {
    // loop over all jets
    for(Candidate *const &jet : *input_collection)
    {
      if(jet->Momentum.Pt() <= fJetPTMin) continue;

      // loop over all constituents
      for(Candidate *const &constituent : jet->GetCandidates())
      {
        // set the IsConstituent flag
        constituent->IsConstituent = 1;
      }
    }
  }

  // loop over all constituent input arrays
  for(const auto &[input_collection, output_collection] : fInputMap)
  {
    // loop over all constituents
    for(Candidate *const &constituent : *input_collection)
    {
      // check the IsConstituent flag
      if(constituent->IsConstituent)
        output_collection->emplace_back(constituent);
    }
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("ConstituentFilter", ConstituentFilter);
