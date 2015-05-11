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

ImpactParameterSmearing::ImpactParameterSmearing() :
  fFormula(0), fItInputArray(0)
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

  fInputArray = ImportArray(GetString("InputArray", "TrackMerger/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void ImpactParameterSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void ImpactParameterSmearing::Process()
{
  Candidate *candidate, *particle, *mother;
  Double_t xd, yd, zd, dxy, sx, sy, sz, ddxy;
  Double_t pt, eta, px, py, phi, e;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {

    // take momentum before smearing (otherwise apply double smearing on dxy)
    particle = static_cast<Candidate*>(candidate->GetCandidates()->At(0));

    const TLorentzVector &candidateMomentum = particle->Momentum;

    eta = candidateMomentum.Eta();
    pt = candidateMomentum.Pt();
    phi = candidateMomentum.Phi();
    e = candidateMomentum.E();
    
    px = candidateMomentum.Px();
    py = candidateMomentum.Py();

    // calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
    xd =  candidate->Xd;
    yd =  candidate->Yd;
    zd =  candidate->Zd;

    // calculate smeared values
    sx = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));
    sy = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));
    sz = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));

    xd += sx;
    yd += sy;
    zd += sz;

    // calculate impact parameter (after-smearing)
    dxy = (xd*py - yd*px)/pt;

    ddxy = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));

    // fill smeared values in candidate
    mother = candidate;

    candidate = static_cast<Candidate*>(candidate->Clone());
    candidate->Xd = xd;
    candidate->Yd = yd;
    candidate->Zd = zd;

    candidate->Dxy = dxy;
    candidate->SDxy = ddxy;

    candidate->AddCandidate(mother);
    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
