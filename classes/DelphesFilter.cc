
/** \class DelphesFilter
 *
 *  Class simplifying classification and subarrays handling
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesFilter.h"

#include "ExRootAnalysis/ExRootClassifier.h"

#include <sstream>
#include <stdexcept>

using namespace std;

typedef map<Int_t, CandidatesCollection> TCategoryMap;
typedef map<ExRootClassifier *, pair<Bool_t, TCategoryMap> > TClassifierMap;

//------------------------------------------------------------------------------

DelphesFilter::DelphesFilter(const CandidatesCollection &collection) : fCollection(collection) {}

//------------------------------------------------------------------------------

void DelphesFilter::Reset(ExRootClassifier *classifier)
{
  TClassifierMap::iterator itMap;
  TCategoryMap::iterator itSubMap;
  if(classifier)
  {
    itMap = fMap.find(classifier);
    if(itMap != fMap.end())
    {
      itMap->second.first = kTRUE;
      for(itSubMap = itMap->second.second.begin();
        itSubMap != itMap->second.second.end(); ++itSubMap)
      {
        itSubMap->second->clear();
      }
    }
  }
  else
  {
    for(itMap = fMap.begin(); itMap != fMap.end(); ++itMap)
    {
      itMap->second.first = kTRUE;
      for(itSubMap = itMap->second.second.begin();
        itSubMap != itMap->second.second.end(); ++itSubMap)
      {
        itSubMap->second->clear();
      }
    }
  }
}

//------------------------------------------------------------------------------

CandidatesCollection DelphesFilter::GetSubArray(ExRootClassifier *classifier, Int_t category)
{
  Int_t result;
  TCategoryMap::iterator itSubMap;
  pair<TCategoryMap::iterator, bool> pairSubMap;
  pair<TClassifierMap::iterator, bool> pairMap;

  TClassifierMap::iterator itMap = fMap.find(classifier);
  if(itMap == fMap.end())
  {
    pairMap = fMap.insert(make_pair(classifier, make_pair(kTRUE, TCategoryMap())));
    if(!pairMap.second)
    {
      throw runtime_error("can't insert category map");
    }

    itMap = pairMap.first;
  }

  if(itMap->second.first)
  {
    itMap->second.first = kFALSE;
    for(const auto &element : *fCollection)
    {
      result = classifier->GetCategory(element);
      if(result < 0) continue;
      itSubMap = itMap->second.second.find(result);
      if(itSubMap == itMap->second.second.end())
      {
        pairSubMap = itMap->second.second.insert(make_pair(result, std::make_shared<std::vector<Candidate *> >()));
        if(!pairSubMap.second)
        {
          throw runtime_error("can't insert category");
        }

        itSubMap = pairSubMap.first;
      }
      itSubMap->second->emplace_back(element);
    }
  }

  itSubMap = itMap->second.second.find(category);
  return (itSubMap != itMap->second.second.end()) ? itSubMap->second : 0;
}

//------------------------------------------------------------------------------
