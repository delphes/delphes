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

#include "modules/ImpactParameterSmearing.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

ImpactParameterSmearing::ImpactParameterSmearing() :
  fFormula(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

ImpactParameterSmearing::~ImpactParameterSmearing()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void ImpactParameterSmearing::Init()
{
  // read resolution formula

  fFormula->Compile(GetString("ResolutionFormula", "0.0"));

  // import input array
  ImportArray(GetString("InputArray", "TrackMerger/tracks"), fInputArray);
  // create output array
  ExportArray(fOutputArray, GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void ImpactParameterSmearing::Finish()
{
}

//------------------------------------------------------------------------------

void ImpactParameterSmearing::Process()
{
  Double_t xd, yd, zd, d0, sx, sy, sz, dd0;
  Double_t pt, eta, px, py, phi, e;

  for(const auto &candidate : *fInputArray) //TODO: ensure const safety
  {

    // take momentum before smearing (otherwise apply double smearing on d0)
    auto *particle = static_cast<Candidate *>(candidate.GetCandidates().at(0));

    const auto &candidateMomentum = particle->Momentum;

    eta = candidateMomentum.Eta();
    pt = candidateMomentum.Pt();
    phi = candidateMomentum.Phi();
    e = candidateMomentum.E();

    px = candidateMomentum.Px();
    py = candidateMomentum.Py();

    // calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
    xd = candidate.Xd;
    yd = candidate.Yd;
    zd = candidate.Zd;

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
    auto new_candidate = candidate;
    new_candidate.Xd = xd;
    new_candidate.Yd = yd;
    new_candidate.Zd = zd;

    new_candidate.D0 = d0;
    new_candidate.ErrorD0 = dd0;

    new_candidate.AddCandidate(const_cast<Candidate *>(&candidate)); // ensure parentage
    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------
