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

/** \class AngularSmearing
 *
 *  Performs transverse angular resolution smearing.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/AngularSmearing.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

AngularSmearing::AngularSmearing() :
  fFormulaEta(0), fFormulaPhi(0), fItInputArray(0)
{
  fFormulaEta = new DelphesFormula;
  fFormulaPhi = new DelphesFormula;
}

//------------------------------------------------------------------------------

AngularSmearing::~AngularSmearing()
{
  if(fFormulaEta) delete fFormulaEta;
  if(fFormulaPhi) delete fFormulaPhi;
}

//------------------------------------------------------------------------------

void AngularSmearing::Init()
{
  // read resolution formula

  fFormulaEta->Compile(GetString("EtaResolutionFormula", "0.0"));
  fFormulaPhi->Compile(GetString("PhiResolutionFormula", "0.0"));

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void AngularSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void AngularSmearing::Process()
{
  Candidate *candidate, *mother;
  Double_t pt, eta, phi, e, m;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = candidateMomentum.Eta();
    phi = candidateMomentum.Phi();
    pt = candidateMomentum.Pt();
    e = candidateMomentum.E();
    m = candidateMomentum.M();

    // apply smearing formula for eta,phi
    eta = gRandom->Gaus(eta, fFormulaEta->Eval(pt, eta, phi, e, candidate));
    phi = gRandom->Gaus(phi, fFormulaPhi->Eval(pt, eta, phi, e, candidate));

    if(pt <= 0.0) continue;

    mother = candidate;
    candidate = static_cast<Candidate *>(candidate->Clone());
    candidate->Momentum.SetPtEtaPhiM(pt, eta, phi, m);
    candidate->AddCandidate(mother);

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
