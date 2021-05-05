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
 *  Performs time smearing.
 *
 *  \author M. Selvaggi - CERN
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
  fItTrackInputArray(0), fResolutionFormula(0)
{
	fResolutionFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

TimeSmearing::~TimeSmearing()
{
	if(fResolutionFormula) delete fResolutionFormula;
}

//------------------------------------------------------------------------------

void TimeSmearing::Init()
{
  // read resolution formula

  // read time resolution formula in seconds
  fResolutionFormula->Compile(GetString("TimeResolution", "30e-12"));

  // import track input array
  fTrackInputArray = ImportArray(GetString("TrackInputArray", "MuonMomentumSmearing/muons"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();


  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void TimeSmearing::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
}

//------------------------------------------------------------------------------

void TimeSmearing::Process()
{
  Candidate *candidate, *mother;
  Double_t tf_smeared, tf;
  Double_t eta, energy;
  Double_t timeResolution;

  const Double_t c_light = 2.99792458E8;

  fItTrackInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItTrackInputArray->Next())))
  {    // converting to meters

    const TLorentzVector &candidateFinalPosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;

    tf = candidateFinalPosition.T() * 1.0E-3 / c_light;
    eta = candidateMomentum.Eta();
    energy = candidateMomentum.E();

    // apply smearing formula
    timeResolution = fResolutionFormula->Eval(0.0, eta, 0.0, energy);
    tf_smeared = gRandom->Gaus(tf, timeResolution);

    mother = candidate;
    candidate = static_cast<Candidate *>(candidate->Clone());

    candidate->Position.SetT(tf_smeared * 1.0E3 * c_light);
    candidate->ErrorT = timeResolution * 1.0E3 * c_light;

    candidate->AddCandidate(mother);
    fOutputArray->Add(candidate);
  }
}
