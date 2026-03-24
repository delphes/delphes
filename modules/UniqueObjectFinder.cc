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

class UniqueObjectFinder: public DelphesModule
{
public:
  explicit UniqueObjectFinder(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fUseUniqueID(Steer<bool>("UseUniqueID", false)) // use GetUniqueID algorithm to find unique objects (faster than the default Overlaps method)
  {
  }

  void Init() override
  {
    // import arrays with output from other modules
    for(const std::pair<std::string, std::string> &collectionsLabel :
      Steer<std::vector<std::pair<std::string, std::string> > >("InputArray"))
      fInputMap.emplace_back(std::make_pair(
        ImportArray(collectionsLabel.first),
        ExportArray(collectionsLabel.second)));
  }
  void Process() override;

private:
  using InputMap = std::vector<std::pair<CandidatesCollection, CandidatesCollection> >;
  bool Unique(Candidate *candidate, InputMap::const_iterator itInputMap);

  const bool fUseUniqueID;

  InputMap fInputMap; //!
};

//------------------------------------------------------------------------------

void UniqueObjectFinder::Process()
{
  for(const auto &[inputCollection, outputCollection] : fInputMap)
    outputCollection->clear();

  for(auto itInputMap = fInputMap.cbegin(); itInputMap != fInputMap.cend(); ++itInputMap) // loop over all input arrays
  {
    for(Candidate *const &candidate : *(itInputMap->first)) // loop over all candidates
    {
      if(Unique(candidate, itInputMap))
        itInputMap->second->emplace_back(candidate);
    }
  }
}

//------------------------------------------------------------------------------

bool UniqueObjectFinder::Unique(Candidate *candidate, InputMap::const_iterator itInputMap)
{
  for(auto previousItInputMap = fInputMap.cbegin(); previousItInputMap != itInputMap; ++previousItInputMap)
  { // loop over previous arrays
    for(Candidate *const &previousCandidate : *(previousItInputMap->second)) // loop over all candidates
      if(fUseUniqueID)
      {
        if(candidate->GetUniqueID() == previousCandidate->GetUniqueID()) return false;
      }
      else if(candidate->Overlaps(previousCandidate))
        return false;
  }
  return true;
}

//------------------------------------------------------------------------------

REGISTER_MODULE("UniqueObjectFinder", UniqueObjectFinder);
