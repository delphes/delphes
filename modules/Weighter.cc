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

class Weighter: public DelphesModule
{
public:
  explicit Weighter(const DelphesParameters &moduleParams) : DelphesModule(moduleParams)
  {
    for(const std::pair<std::vector<int>, double> &weightValue :
      Steer<std::vector<std::pair<std::vector<int>, double> > >("Weight"))
    {
      if(const size_t sizeCodes = weightValue.first.size(); sizeCodes < 1 || sizeCodes > 4)
        throw std::runtime_error("only 1, 2, 3 or 4 PDG codes can be specified per weight");
      TIndexStruct index;
      std::memset(index.codes, 0, 4 * sizeof(int));
      size_t j = 0;
      for(const int &code : weightValue.first)
      {
        index.codes[j++] = code;
        fWeightSet.insert(code);
      }
      std::sort(index.codes, index.codes + 4);
      fWeightMap[index] = weightValue.second;
    }

    TIndexStruct index;
    std::memset(index.codes, 0, 4 * sizeof(int));
    fWeightMap[index] = 1.0;
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "Delphes/allParticles"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "weight"));
  }
  void Process() override;

private:
#if !defined(__CINT__) && !defined(__CLING__)
  struct TIndexStruct
  {
    int codes[4];
    bool operator<(const TIndexStruct &value) const
    {
      for(size_t i = 0; i < 4; ++i)
        if(codes[i] != value.codes[i]) return codes[i] < value.codes[i];
      return false;
    }
  };

  std::set<int> fWeightSet, fCodeSet;
  std::map<TIndexStruct, double> fWeightMap;
#endif

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void Weighter::Process()
{
  fOutputArray->clear();

  DelphesFactory *factory = GetFactory();

  // loop over all particles
  fCodeSet.clear();
  for(Candidate *const &candidate : *fInputArray)
  {
    if(candidate->Status != 3) continue;
    if(fWeightSet.find(candidate->PID) == fWeightSet.end()) continue;

    fCodeSet.insert(candidate->PID);
  }

  // find default weight value
  TIndexStruct index;
  std::memset(index.codes, 0, 4 * sizeof(int));
  double weight = fWeightMap.at(index);

  if(fCodeSet.size() <= 4)
  {
    size_t i = 0;
    for(const int &codeSet : fCodeSet)
      index.codes[i++] = codeSet;

    std::sort(index.codes, index.codes + 4);
    if(fWeightMap.count(index) > 0)
      weight = fWeightMap.at(index);
  }

  Candidate *candidate = factory->NewCandidate();
  candidate->Momentum.SetPtEtaPhiE(weight, 0.0, 0.0, weight);
  fOutputArray->emplace_back(candidate);
}

//------------------------------------------------------------------------------

REGISTER_MODULE("Weighter", Weighter);
