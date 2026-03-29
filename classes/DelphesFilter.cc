
/** \class DelphesFilter
 *
 *  Class simplifying classification and subarrays handling
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesFilter.h"

#include "ExRootAnalysis/ExRootClassifier.h"

#include <stdexcept>

DelphesFilter::DelphesFilter(const CandidatesCollection &collection) : fCollection(collection) {}

//------------------------------------------------------------------------------

void DelphesFilter::Reset(ExRootClassifier *classifier)
{
  if(classifier)
  {
    TClassifierMap::iterator itMap = fMap.find(classifier);
    if(itMap != fMap.end())
    {
      itMap->second.first = true;
      for(TCategoryMap::iterator itSubMap = itMap->second.second.begin(); itSubMap != itMap->second.second.end(); ++itSubMap)
        itSubMap->second.clear();
    }
  }
  else
  {
    for(TClassifierMap::iterator itMap = fMap.begin(); itMap != fMap.end(); ++itMap)
    {
      itMap->second.first = true;
      for(TCategoryMap::iterator itSubMap = itMap->second.second.begin(); itSubMap != itMap->second.second.end(); ++itSubMap)
        itSubMap->second.clear();
    }
  }
}

//------------------------------------------------------------------------------

std::vector<Candidate *> DelphesFilter::GetSubArray(ExRootClassifier *classifier, Int_t category)
{
  TClassifierMap::iterator itMap = fMap.find(classifier);
  if(itMap == fMap.end())
  {
    std::pair<TClassifierMap::iterator, bool> pairMap = fMap.insert(
      std::make_pair(classifier, std::make_pair(true, TCategoryMap())));
    if(!pairMap.second)
      throw std::runtime_error("can't insert category map");
    itMap = pairMap.first;
  }

  if(itMap->second.first)
  {
    itMap->second.first = false;
    for(Candidate *const &element : *fCollection)
    {
      int result = classifier->GetCategory(element);
      if(result < 0) continue;
      TCategoryMap::iterator itSubMap = itMap->second.second.find(result);
      if(itSubMap == itMap->second.second.end())
      {
        std::pair<TCategoryMap::iterator, bool> pairSubMap = itMap->second.second.insert(
          std::make_pair(result, std::vector<Candidate *>()));
        if(!pairSubMap.second)
          throw std::runtime_error("can't insert category");

        itSubMap = pairSubMap.first;
      }
      itSubMap->second.emplace_back(element);
    }
  }
  TCategoryMap::iterator itSubMap = itMap->second.second.find(category);
  return (itSubMap != itMap->second.second.end()) ? itSubMap->second : std::vector<Candidate *>{};
}

//------------------------------------------------------------------------------
