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

/** \class PdgCodeFilter
 *
 *  Removes particles with specific PDG codes
 *
 *  \author M. Selvaggi
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"
#include "classes/DelphesModuleFactory.h"

#include <TLorentzVector.h>

using namespace std;

class PdgCodeFilter: public DelphesModule
{
public:
  PdgCodeFilter() = default;

  void Init() override;
  void Process() override;

private:
  Double_t fPTMin; //!
  Bool_t fInvert; //!
  Bool_t fRequireStatus; //!
  Int_t fStatus; //!
  Bool_t fRequireCharge; //!
  Int_t fCharge; //!
  Bool_t fRequireNotPileup; //!

  std::vector<Int_t> fPdgCodes;

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void PdgCodeFilter::Init()
{
  ExRootConfParam param;
  Size_t i, size;

  // PT threshold
  fPTMin = GetDouble("PTMin", 0.0);

  fInvert = GetBool("Invert", false);

  // no pileup
  fRequireNotPileup = GetBool("RequireNotPileup", false);

  fRequireStatus = GetBool("RequireStatus", false);
  fStatus = GetInt("Status", 1);

  fRequireCharge = GetBool("RequireCharge", false);
  fCharge = GetInt("Charge", 1);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));

  param = GetParam("PdgCode");
  size = param.GetSize();

  // read PdgCodes to be filtered out from the data card

  fPdgCodes.clear();
  for(i = 0; i < size; ++i)
  {
    fPdgCodes.push_back(param[i].GetInt());
  }

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void PdgCodeFilter::Process()
{
  fOutputArray->clear();

  Int_t pdgCode;
  Bool_t pass;
  Double_t pt;

  for(const auto &candidate : *fInputArray)
  {
    pdgCode = candidate->PID;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    pt = candidateMomentum.Pt();

    if(pt < fPTMin) continue;
    if(fRequireStatus && (candidate->Status != fStatus)) continue;
    if(fRequireCharge && (candidate->Charge != fCharge)) continue;
    if(fRequireNotPileup && (candidate->IsPU > 0)) continue;

    pass = kTRUE;
    if(find(fPdgCodes.begin(), fPdgCodes.end(), pdgCode) != fPdgCodes.end()) pass = kFALSE;

    if(fInvert) pass = !pass;
    if(pass) fOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("PdgCodeFilter", PdgCodeFilter);
