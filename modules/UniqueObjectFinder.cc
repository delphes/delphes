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
#include "classes/DelphesModuleFactory.h"

#include <TObjArray.h>

#include <vector>

using namespace std;

class UniqueObjectFinder: public DelphesModule
{
public:
  UniqueObjectFinder() = default;

  void Init() override;
  void Process() override;
  void Finish() override;

private:
  Bool_t fUseUniqueID;

  using InputMap = std::vector<std::pair<std::unique_ptr<TIterator>, TObjArray *> >;

  Bool_t Unique(Candidate *candidate, InputMap::iterator itInputMap);

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
  const TObjArray *array = nullptr;

  fInputMap.clear();

  size = param.GetSize();
  for(i = 0; i < size / 2; ++i)
  {
    array = ImportArray(param[i * 2].GetString());
    auto iterator = std::unique_ptr<TIterator>(array->MakeIterator());

    fInputMap.push_back(make_pair(std::move(iterator), ExportArray(param[i * 2 + 1].GetString())));
  }
}

//------------------------------------------------------------------------------

void UniqueObjectFinder::Finish()
{
  fInputMap.clear();
}

//------------------------------------------------------------------------------

void UniqueObjectFinder::Process()
{
  Candidate *candidate = nullptr;
  InputMap::iterator itInputMap;
  TObjArray *array = nullptr;

  // loop over all input arrays
  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    auto &iterator = itInputMap->first;
    array = itInputMap->second;

    // loop over all candidates
    iterator->Reset();
    while((candidate = static_cast<Candidate *>(iterator->Next())))
    {
      if(Unique(candidate, itInputMap))
      {
        array->Add(candidate);
      }
    }
  }
}

//------------------------------------------------------------------------------

Bool_t UniqueObjectFinder::Unique(Candidate *candidate, InputMap::iterator itInputMap)
{
  Candidate *previousCandidate = nullptr;
  InputMap::iterator previousItInputMap;
  TObjArray *array = nullptr;

  // loop over previous arrays
  for(previousItInputMap = fInputMap.begin(); previousItInputMap != itInputMap; ++previousItInputMap)
  {
    array = previousItInputMap->second;
    TIter iterator(array);

    // loop over all candidates
    iterator.Reset();
    while((previousCandidate = static_cast<Candidate *>(iterator.Next())))
    {
      if(fUseUniqueID)
      {
        if(candidate->GetUniqueID() == previousCandidate->GetUniqueID())
        {
          return kFALSE;
        }
      }
      else
      {
        if(candidate->Overlaps(previousCandidate))
        {
          return kFALSE;
        }
      }
    }
  }

  return kTRUE;
}

//------------------------------------------------------------------------------

REGISTER_MODULE("UniqueObjectFinder", UniqueObjectFinder);
