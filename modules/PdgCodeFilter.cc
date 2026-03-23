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

#include <TLorentzVector.h>

using namespace std;

class PdgCodeFilter: public DelphesModule
{
public:
  explicit PdgCodeFilter(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fPTMin(Steer<double>("PTMin", 0.0)), // PT threshold
    fInvert(Steer<bool>("Invert", false)),
    fRequireNotPileup(Steer<bool>("RequireNotPileup", false)), // no pileup
    fRequireStatus(Steer<bool>("RequireStatus", false)),
    fStatus(Steer<int>("Status", 1)),
    fRequireCharge(Steer<bool>("RequireCharge", false)),
    fCharge(Steer<int>("Charge", 1)),
    fPdgCodes(Steer<std::vector<int> >("PdgCode")) // PdgCodes to be filtered out from the data card
  {
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "Delphes/allParticles"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "filteredParticles"));
  }
  void Process() override;

private:
  const double fPTMin; //!
  const bool fInvert; //!
  const bool fRequireNotPileup; //!
  const bool fRequireStatus; //!
  const int fStatus; //!
  const bool fRequireCharge; //!
  const int fCharge; //!

  const std::vector<int> fPdgCodes;

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void PdgCodeFilter::Process()
{
  fOutputArray->clear();

  for(Candidate *const &candidate : *fInputArray)
  {
    const int pdgCode = candidate->PID;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    const double pt = candidateMomentum.Pt();

    if(pt < fPTMin) continue;
    if(fRequireStatus && (candidate->Status != fStatus)) continue;
    if(fRequireCharge && (candidate->Charge != fCharge)) continue;
    if(fRequireNotPileup && (candidate->IsPU > 0)) continue;

    bool pass = std::find(fPdgCodes.begin(), fPdgCodes.end(), pdgCode) == fPdgCodes.end();
    if(fInvert) pass = !pass;
    if(pass) fOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("PdgCodeFilter", PdgCodeFilter);
