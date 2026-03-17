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
#include "classes/DelphesModuleFactory.h"

#include <TObjArray.h>

using namespace std;

class ConstituentFilter: public DelphesModule
{
public:
  ConstituentFilter() = default;

  void Init() override;
  void Process() override;
  void Finish() override;

private:
  Double_t fJetPTMin;

  std::vector<std::unique_ptr<TIterator> > fInputList; //!

  std::map<std::unique_ptr<TIterator>, TObjArray *> fInputMap; //!

  TObjArray *fOutputArray{nullptr}; //!
};

//------------------------------------------------------------------------------

void ConstituentFilter::Init()
{
  ExRootConfParam param;
  Long_t i, size;
  const TObjArray *array;

  fJetPTMin = GetDouble("JetPTMin", 0.0);

  // import input array(s)

  param = GetParam("JetInputArray");
  size = param.GetSize();
  for(i = 0; i < size; ++i)
  {
    array = ImportArray(param[i].GetString());
    fInputList.push_back(std::unique_ptr<TIterator>(array->MakeIterator()));
  }

  param = GetParam("ConstituentInputArray");
  size = param.GetSize();
  for(i = 0; i < size / 2; ++i)
  {
    array = ImportArray(param[i * 2].GetString());
    fInputMap[std::unique_ptr<TIterator>(array->MakeIterator())] = ExportArray(param[i * 2 + 1].GetString());
  }
}

//------------------------------------------------------------------------------

void ConstituentFilter::Finish()
{
  fInputList.clear();
  fInputMap.clear();
}

//------------------------------------------------------------------------------

void ConstituentFilter::Process()
{
  Candidate *jet, *constituent;
  map<std::unique_ptr<TIterator>, TObjArray *>::iterator itInputMap;
  vector<std::unique_ptr<TIterator> >::iterator itInputList;
  TObjArray *array;

  // loop over all jet input arrays
  for(itInputList = fInputList.begin(); itInputList != fInputList.end(); ++itInputList)
  {
    auto &iterator = *itInputList;

    // loop over all jets
    iterator->Reset();
    while((jet = static_cast<Candidate *>(iterator->Next())))
    {
      TIter itConstituents(jet->GetCandidates());

      if(jet->Momentum.Pt() <= fJetPTMin) continue;

      // loop over all constituents
      while((constituent = static_cast<Candidate *>(itConstituents.Next())))
      {
        // set the IsConstituent flag
        constituent->IsConstituent = 1;
      }
    }
  }

  // loop over all constituent input arrays
  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    auto &iterator = itInputMap->first;
    array = itInputMap->second;

    // loop over all constituents
    iterator->Reset();
    while((constituent = static_cast<Candidate *>(iterator->Next())))
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

REGISTER_MODULE("ConstituentFilter", ConstituentFilter);
