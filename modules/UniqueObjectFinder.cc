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

/** \class UniqueObjectFinder
 *
 *  Finds uniquely identified photons, electrons, taus and jets.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"

#include <vector>

using namespace std;

class UniqueObjectFinder: public DelphesModule
{
public:
  UniqueObjectFinder() = default;

  void Init() override;
  void Process() override;

private:
  Bool_t fUseUniqueID;

  using InputMap = std::vector<std::pair<CandidatesCollection, CandidatesCollection> >;

  Bool_t Unique(Candidate *candidate, InputMap::const_iterator itInputMap);

  InputMap fInputMap; //!
};

//------------------------------------------------------------------------------

void UniqueObjectFinder::Init()
{
  // use GetUniqueID algorithm to find unique objects (faster than the default Overlaps method)
  fUseUniqueID = GetBool("UseUniqueID", false);

  // import arrays with output from other modules

  ExRootConfParam param = GetParam("InputArray");
  Long_t i, size;

  fInputMap.clear();

  size = param.GetSize();
  for(i = 0; i < size / 2; ++i)
    fInputMap.emplace_back(std::make_pair(ImportArray(param[i * 2].GetString()), ExportArray(param[i * 2 + 1].GetString())));
}

//------------------------------------------------------------------------------

void UniqueObjectFinder::Process()
{
  for(const auto &[input_collection, output_collection] : fInputMap)
    output_collection->clear();

  for(auto itInputMap = fInputMap.cbegin(); itInputMap != fInputMap.cend(); ++itInputMap) // loop over all input arrays
  {
    for(const auto &candidate : *(itInputMap->first)) // loop over all candidates
    {
      if(Unique(candidate, itInputMap))
        itInputMap->second->emplace_back(candidate);
    }
  }
}

//------------------------------------------------------------------------------

Bool_t UniqueObjectFinder::Unique(Candidate *candidate, InputMap::const_iterator itInputMap)
{
  for(auto previousItInputMap = fInputMap.cbegin(); previousItInputMap != itInputMap; ++previousItInputMap)
  { // loop over previous arrays
    for(const auto &previousCandidate : *(previousItInputMap->second)) // loop over all candidates
      if(fUseUniqueID)
      {
        if(candidate->GetUniqueID() == previousCandidate->GetUniqueID()) return kFALSE;
      }
      else
      {
        if(candidate->Overlaps(previousCandidate)) return kFALSE;
      }
  }
  return kTRUE;
}

//------------------------------------------------------------------------------

REGISTER_MODULE("UniqueObjectFinder", UniqueObjectFinder);
