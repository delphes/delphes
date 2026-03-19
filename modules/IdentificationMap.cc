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

/** \class IdentificationMap
 *
 *  Converts particles with some PDG code into another particle,
 *  according to parametrized probability.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

using namespace std;

class IdentificationMap: public DelphesModule
{
public:
  IdentificationMap() = default;

  void Init() override;
  void Process() override;
  void Finish() override;

private:
  typedef std::multimap<Int_t, std::pair<Int_t, std::unique_ptr<DelphesFormula> > > TMisIDMap; //!

  TMisIDMap fEfficiencyMap; //!

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void IdentificationMap::Init()
{
  TMisIDMap::iterator itEfficiencyMap;
  ExRootConfParam param;
  Int_t i, size, pdg;

  // read efficiency formulas
  param = GetParam("EfficiencyFormula");
  size = param.GetSize();

  fEfficiencyMap.clear();
  for(i = 0; i < size / 3; ++i)
  {
    std::unique_ptr<DelphesFormula> formula = std::make_unique<DelphesFormula>();
    formula->Compile(param[i * 3 + 2].GetString());
    pdg = param[i * 3].GetInt();
    fEfficiencyMap.insert(make_pair(pdg, make_pair(param[i * 3 + 1].GetInt(), std::move(formula))));
  }

  // set default efficiency formula
  itEfficiencyMap = fEfficiencyMap.find(0);
  if(itEfficiencyMap == fEfficiencyMap.end())
  {
    std::unique_ptr<DelphesFormula> formula = std::make_unique<DelphesFormula>();
    formula->Compile("1.0");

    fEfficiencyMap.insert(make_pair(0, make_pair(0, std::move(formula))));
  }

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void IdentificationMap::Finish()
{
  fEfficiencyMap.clear();
}

//------------------------------------------------------------------------------

void IdentificationMap::Process()
{
  fOutputArray->clear();

  Double_t pt, eta, phi, e;
  TMisIDMap::iterator itEfficiencyMap;
  pair<TMisIDMap::iterator, TMisIDMap::iterator> range;
  Int_t pdgCodeIn, pdgCodeOut, charge;

  Double_t p, r, total;

  for(Candidate *const &candidate : *fInputArray)
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();
    pt = candidateMomentum.Pt();
    e = candidateMomentum.E();

    pdgCodeIn = candidate->PID;
    charge = candidate->Charge;

    // first check that PID of this particle is specified in the map
    // otherwise, look for PID = 0

    itEfficiencyMap = fEfficiencyMap.find(pdgCodeIn);

    range = fEfficiencyMap.equal_range(pdgCodeIn);
    if(range.first == range.second) range = fEfficiencyMap.equal_range(-pdgCodeIn);
    if(range.first == range.second) range = fEfficiencyMap.equal_range(0);

    r = gRandom->Uniform();
    total = 0.0;

    // loop over sub-map for this PID
    for(TMisIDMap::iterator it = range.first; it != range.second; ++it)
    {
      std::unique_ptr<DelphesFormula> &formula = (it->second).second;
      pdgCodeOut = (it->second).first;

      p = formula->Eval(pt, eta, phi, e);

      if(total <= r && r < total + p)
      {
        // change PID of particle
        Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());
        if(pdgCodeOut != 0) new_candidate->PID = charge * pdgCodeOut;
        fOutputArray->emplace_back(new_candidate);
        break;
      }

      total += p;
    }
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("IdentificationMap", IdentificationMap);
