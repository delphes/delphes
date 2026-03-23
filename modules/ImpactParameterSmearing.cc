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

/** \class ImpactParameterSmearing
 *
 *  Performs transverse impact parameter smearing.
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

class ImpactParameterSmearing: public DelphesModule
{
public:
  explicit ImpactParameterSmearing(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fFormula(std::make_unique<DelphesFormula>())
  {
    // read resolution formula
    fFormula->Compile(Steer<std::string>("ResolutionFormula", "0.0"));
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "TrackMerger/tracks"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "tracks"));
  }
  void Process() override;

private:
  const std::unique_ptr<DelphesFormula> fFormula; //!

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void ImpactParameterSmearing::Process()
{
  fOutputArray->clear();

  for(Candidate *const &candidate : *fInputArray)
  {
    // take momentum before smearing (otherwise apply double smearing on d0)
    Candidate *particle = static_cast<Candidate *>(candidate->GetCandidates().at(0));

    const TLorentzVector &candidateMomentum = particle->Momentum;

    const double pt = candidateMomentum.Pt();
    const double eta = candidateMomentum.Eta();
    const double phi = candidateMomentum.Phi();
    const double e = candidateMomentum.E();

    const double px = candidateMomentum.Px();
    const double py = candidateMomentum.Py();

    // calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
    double xd = candidate->Xd, yd = candidate->Yd, zd = candidate->Zd;

    // calculate smeared values
    const double sx = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));
    const double sy = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));
    const double sz = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));

    xd += sx;
    yd += sy;
    zd += sz;

    // calculate impact parameter (after-smearing)
    const double d0 = (xd * py - yd * px) / pt;
    const double dd0 = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));

    // fill smeared values in candidate
    Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());
    new_candidate->Xd = xd;
    new_candidate->Yd = yd;
    new_candidate->Zd = zd;

    new_candidate->D0 = d0;
    new_candidate->ErrorD0 = dd0;

    new_candidate->AddCandidate(candidate);
    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("ImpactParameterSmearing", ImpactParameterSmearing);
