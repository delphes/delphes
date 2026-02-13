
/** \class ExRootSTLVectorFilter
 *
 *  Class simplifying classification and subarrays handling
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *          L. Forthomme - AGH, Krakow
 *
 */

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootSTLVectorFilter.h"

#include "classes/DelphesClasses.h" //TODO: should not depend on Delphes objects... make this templated?

#include "TObjArray.h"
#include "TSeqCollection.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

ExRootSTLVectorFilter::ExRootSTLVectorFilter(const std::vector<Candidate> &collection) :
  fCollection(collection)
{
}

//------------------------------------------------------------------------------

void ExRootSTLVectorFilter::Reset(ExRootClassifier *classifier)
{
  if(classifier)
  {
    if(auto itMap = fMap.find(classifier); itMap != fMap.end())
    {
      itMap->second.first = kTRUE;
      for(auto itSubMap = itMap->second.second.begin(); itSubMap != itMap->second.second.end(); ++itSubMap)
        itSubMap->second.clear();
    }
  }
  else
  {
    for(auto itMap = fMap.begin(); itMap != fMap.end(); ++itMap)
    {
      itMap->second.first = kTRUE;
      for(auto itSubMap = itMap->second.second.begin(); itSubMap != itMap->second.second.end(); ++itSubMap)
        itSubMap->second.clear();
    }
  }
}

//------------------------------------------------------------------------------

std::vector<Candidate> ExRootSTLVectorFilter::GetSubArray(ExRootClassifier *classifier, Int_t category)
{
  auto itMap = fMap.find(classifier);
  if(itMap == fMap.end()) // classifier was not found
  {
    auto pairMap = fMap.insert(make_pair(classifier, make_pair(kTRUE, TCategoryMap())));
    if(!pairMap.second)
      throw runtime_error("can't insert category map");
    itMap = pairMap.first;
  }

  if(itMap->second.first)
  {
    itMap->second.first = kFALSE;
    for(const auto &element : fCollection)
    {
      const auto result = classifier->GetCategory(const_cast<TObject *>(static_cast<const TObject *>(&element)));
      if(result < 0) continue;
      auto itSubMap = itMap->second.second.find(result);
      if(itSubMap == itMap->second.second.end())
      {
        auto pairSubMap = itMap->second.second.insert(make_pair(result, std::vector<Candidate>(fCollection.size())));
        if(!pairSubMap.second)
        {
          throw runtime_error("can't insert category");
        }

        itSubMap = pairSubMap.first;
      }
      itSubMap->second.emplace_back(element);
    }
  }

  auto itSubMap = itMap->second.second.find(category);
  return (itSubMap != itMap->second.second.end()) ? itSubMap->second : std::vector<Candidate>{};
}

//------------------------------------------------------------------------------
