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

/** \class TimeSmearing
 *
 *  Performs transverse momentum resolution smearing.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/TimeSmearing.h"

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

TimeSmearing::TimeSmearing() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

TimeSmearing::~TimeSmearing()
{
}

//------------------------------------------------------------------------------

void TimeSmearing::Init()
{
  // read resolution formula

  fTimeResolution = GetDouble("TimeResolution", 1.0E-10);
  // import input array

  fInputArray = ImportArray(GetString("InputArray", "MuonMomentumSmearing/muons"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "muons"));
}

//------------------------------------------------------------------------------

void TimeSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void TimeSmearing::Process()
{
  Candidate *candidate, *mother;
  Double_t ti, tf_smeared, tf;
  const Double_t c_light = 2.99792458E8;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    const TLorentzVector &candidateInitialPosition = candidate->InitialPosition;
    const TLorentzVector &candidateFinalPosition = candidate->Position;

    ti = candidateInitialPosition.T() * 1.0E-3 / c_light;
    tf = candidateFinalPosition.T() * 1.0E-3 / c_light;

    // apply smearing formula
    tf_smeared = gRandom->Gaus(tf, fTimeResolution);
    ti = ti + tf_smeared - tf;
    tf = tf_smeared;

    mother = candidate;
    candidate = static_cast<Candidate *>(candidate->Clone());
    candidate->InitialPosition.SetT(ti * 1.0E3 * c_light);
    candidate->Position.SetT(tf * 1.0E3 * c_light);

    candidate->ErrorT = fTimeResolution * 1.0E3 * c_light;

    candidate->AddCandidate(mother);

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
