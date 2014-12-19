
/** \class ExRootFilter
 *
 *  Class simplifying classification and subarrays handling
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TSeqCollection.h"
#include "TObjArray.h"

#include <iostream>
#include <stdexcept>
#include <sstream>

using namespace std;

typedef map<Int_t, TObjArray*> TCategoryMap;
typedef map<ExRootClassifier*, pair<Bool_t, TCategoryMap> > TClassifierMap;

ExRootFilter::ExRootFilter(const TSeqCollection *collection) :
  fCollection(collection)
{
  fIter = fCollection->MakeIterator();
}

//------------------------------------------------------------------------------

ExRootFilter::~ExRootFilter()
{
  TClassifierMap::iterator itMap;
  TCategoryMap::iterator itSubMap;
  for(itMap = fMap.begin(); itMap != fMap.end(); ++itMap)
  {
    for(itSubMap = itMap->second.second.begin();
        itSubMap != itMap->second.second.end(); ++itSubMap)
    {
      delete (itSubMap->second);
    }
  }

  delete fIter;
}

//------------------------------------------------------------------------------

void ExRootFilter::Reset(ExRootClassifier *classifier)
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
        itSubMap->second->Clear();
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
        itSubMap->second->Clear();
      }
    }
  }
}

//------------------------------------------------------------------------------

TObjArray *ExRootFilter::GetSubArray(ExRootClassifier *classifier, Int_t category)
{
  Int_t result;
  TObject *element;
  TObjArray *array;
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
    fIter->Reset();
    while((element = fIter->Next()) != 0)
    {
      result = classifier->GetCategory(element);
      if(result < 0) continue;
      itSubMap = itMap->second.second.find(result);
      if(itSubMap == itMap->second.second.end())
      {
        array = new TObjArray(fCollection->GetSize());
        pairSubMap = itMap->second.second.insert(make_pair(result, array));
        if(!pairSubMap.second)
        {
          throw runtime_error("can't insert category");
        }

        itSubMap = pairSubMap.first;
      }
      itSubMap->second->Add(element);
    }
  }

  itSubMap = itMap->second.second.find(category);
  return (itSubMap != itMap->second.second.end()) ? itSubMap->second : 0;
}

//------------------------------------------------------------------------------

