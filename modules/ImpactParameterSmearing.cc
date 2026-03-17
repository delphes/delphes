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
#include "classes/DelphesModuleFactory.h"

#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TRandom3.h>

using namespace std;

class ImpactParameterSmearing: public DelphesModule
{
public:
  ImpactParameterSmearing() : fFormula(std::make_unique<DelphesFormula>()) {}

  void Init() override;
  void Process() override;

private:
  const std::unique_ptr<DelphesFormula> fFormula; //!

  const TObjArray *fInputArray{nullptr}; //!
  std::unique_ptr<TIterator> fItInputArray; //!

  TObjArray *fOutputArray{nullptr}; //!
};

//------------------------------------------------------------------------------

void ImpactParameterSmearing::Init()
{
  // read resolution formula
  fFormula->Compile(GetString("ResolutionFormula", "0.0"));

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "TrackMerger/tracks"));
  fItInputArray.reset(fInputArray->MakeIterator());

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void ImpactParameterSmearing::Process()
{
  Candidate *candidate = nullptr, *particle = nullptr, *mother = nullptr;
  Double_t xd, yd, zd, d0, sx, sy, sz, dd0;
  Double_t pt, eta, px, py, phi, e;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {

    // take momentum before smearing (otherwise apply double smearing on d0)
    particle = static_cast<Candidate *>(candidate->GetCandidates()->At(0));

    const TLorentzVector &candidateMomentum = particle->Momentum;

    eta = candidateMomentum.Eta();
    pt = candidateMomentum.Pt();
    phi = candidateMomentum.Phi();
    e = candidateMomentum.E();

    px = candidateMomentum.Px();
    py = candidateMomentum.Py();

    // calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
    xd = candidate->Xd;
    yd = candidate->Yd;
    zd = candidate->Zd;

    // calculate smeared values
    sx = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));
    sy = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));
    sz = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));

    xd += sx;
    yd += sy;
    zd += sz;

    // calculate impact parameter (after-smearing)
    d0 = (xd * py - yd * px) / pt;

    dd0 = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));

    // fill smeared values in candidate
    mother = candidate;

    candidate = static_cast<Candidate *>(candidate->Clone());
    candidate->Xd = xd;
    candidate->Yd = yd;
    candidate->Zd = zd;

    candidate->D0 = d0;
    candidate->ErrorD0 = dd0;

    candidate->AddCandidate(mother);
    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("ImpactParameterSmearing", ImpactParameterSmearing);
