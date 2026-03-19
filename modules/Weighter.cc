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

/** \class Weighter
 *
 *  Apply a weight depending on PDG code.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesModule.h"

#include <set>

using namespace std;

class Weighter: public DelphesModule
{
public:
  Weighter() = default;

  void Init() override;
  void Process() override;

private:
#if !defined(__CINT__) && !defined(__CLING__)
  struct TIndexStruct
  {
    Int_t codes[4];
    bool operator<(const TIndexStruct &value) const
    {
      for(size_t i = 0; i < 4; ++i)
        if(codes[i] != value.codes[i]) return codes[i] < value.codes[i];
      return false;
    }
  };

  std::set<Int_t> fWeightSet, fCodeSet;
  std::map<TIndexStruct, Double_t> fWeightMap;
#endif

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void Weighter::Init()
{
  ExRootConfParam param, paramCodes;
  Int_t i, j, size, sizeCodes;
  Int_t code;
  TIndexStruct index;
  Double_t weight;

  fWeightSet.clear();

  // set default weight value
  fWeightMap.clear();
  memset(index.codes, 0, 4 * sizeof(Int_t));
  fWeightMap[index] = 1.0;

  // read weights
  param = GetParam("Weight");
  size = param.GetSize();
  for(i = 0; i < size / 2; ++i)
  {
    paramCodes = param[i * 2];
    sizeCodes = paramCodes.GetSize();
    weight = param[i * 2 + 1].GetDouble();

    if(sizeCodes < 1 || sizeCodes > 4)
    {
      throw runtime_error("only 1, 2, 3 or 4 PDG codes can be specified per weight");
    }

    memset(index.codes, 0, 4 * sizeof(Int_t));

    for(j = 0; j < sizeCodes; ++j)
    {
      code = paramCodes[j].GetInt();
      index.codes[j] = code;
      fWeightSet.insert(code);
    }

    sort(index.codes, index.codes + 4);

    fWeightMap[index] = weight;
  }

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "weight"));
}

//------------------------------------------------------------------------------

void Weighter::Process()
{
  fOutputArray->clear();

  Int_t i;
  TIndexStruct index;
  Double_t weight;
  set<Int_t>::iterator itCodeSet;
  map<TIndexStruct, Double_t>::iterator itWeightMap;

  DelphesFactory *factory = GetFactory();

  // loop over all particles
  fCodeSet.clear();
  for(const auto &candidate : *fInputArray)
  {
    if(candidate->Status != 3) continue;

    if(fWeightSet.find(candidate->PID) == fWeightSet.end()) continue;

    fCodeSet.insert(candidate->PID);
  }

  // find default weight value
  memset(index.codes, 0, 4 * sizeof(Int_t));
  itWeightMap = fWeightMap.find(index);
  weight = itWeightMap->second;

  if(fCodeSet.size() <= 4)
  {
    i = 0;
    for(itCodeSet = fCodeSet.begin(); itCodeSet != fCodeSet.end(); ++itCodeSet)
    {
      index.codes[i] = *itCodeSet;
      ++i;
    }

    sort(index.codes, index.codes + 4);

    itWeightMap = fWeightMap.find(index);
    if(itWeightMap != fWeightMap.end())
    {
      weight = itWeightMap->second;
    }
  }

  auto *candidate = factory->NewCandidate();
  candidate->Momentum.SetPtEtaPhiE(weight, 0.0, 0.0, weight);
  fOutputArray->emplace_back(candidate);
}

//------------------------------------------------------------------------------

REGISTER_MODULE("Weighter", Weighter);
