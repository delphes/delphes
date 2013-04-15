
/** \class ConstituentFilter
 *
 *  Drops all input objects that are not constituents of any jet.
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/ConstituentFilter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

ConstituentFilter::ConstituentFilter()
{
}

//------------------------------------------------------------------------------

ConstituentFilter::~ConstituentFilter()
{
}

//------------------------------------------------------------------------------

void ConstituentFilter::Init()
{
  ExRootConfParam param;
  Long_t i, size;
  const TObjArray *array;
  TIterator *iterator;

  fJetPTMin = GetDouble("JetPTMin", 0.0);

  // import input array(s)

  param = GetParam("JetInputArray");
  size = param.GetSize();
  for(i = 0; i < size; ++i)
  {
    array = ImportArray(param[i].GetString());
    iterator = array->MakeIterator();

    fInputList.push_back(iterator);
  }

  param = GetParam("ConstituentInputArray");
  size = param.GetSize();
  for(i = 0; i < size/2; ++i)
  {
    array = ImportArray(param[i*2].GetString());
    iterator = array->MakeIterator();

    fInputMap[iterator] = ExportArray(param[i*2 + 1].GetString());;
  }
}

//------------------------------------------------------------------------------

void ConstituentFilter::Finish()
{
  map< TIterator *, TObjArray * >::iterator itInputMap;
  vector< TIterator * >::iterator itInputList;
  TIterator *iterator;

  for(itInputList = fInputList.begin(); itInputList != fInputList.end(); ++itInputList)
  {
    iterator = *itInputList;
    if(iterator) delete iterator;
  }

  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    iterator = itInputMap->first;
    if(iterator) delete iterator;
  }
}

//------------------------------------------------------------------------------

void ConstituentFilter::Process()
{
  Candidate *jet, *constituent;
  map< TIterator *, TObjArray * >::iterator itInputMap;
  vector< TIterator * >::iterator itInputList;
  TIterator *iterator;
  TObjArray *array;

  // loop over all jet input arrays
  for(itInputList = fInputList.begin(); itInputList != fInputList.end(); ++itInputList)
  {
    iterator = *itInputList;

    // loop over all jets
    iterator->Reset();
    while((jet = static_cast<Candidate*>(iterator->Next())))
    {
      TIter itConstituents(jet->GetCandidates());

      if(jet->Momentum.Pt() <= fJetPTMin) continue;

      // loop over all constituents
      while((constituent = static_cast<Candidate*>(itConstituents.Next())))
      {
        // set the IsConstituent flag
        constituent->IsConstituent = 1;
      }
    }
  }

  // loop over all constituent input arrays
  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    iterator = itInputMap->first;
    array = itInputMap->second;

    // loop over all constituents
    iterator->Reset();
    while((constituent = static_cast<Candidate*>(iterator->Next())))
    {
      // check the IsConstituent flag
      if(constituent->IsConstituent)
      {
        array->Add(constituent);
      }
    }
  }
}

//------------------------------------------------------------------------------
