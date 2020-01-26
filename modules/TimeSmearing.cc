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
  fItInputArray(0), fResolutionFormula(0)
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
  // read time resolution formula in seconds
  fResolutionFormula->Compile(GetString("TimeResolution", "30e-12"));

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "TrackMerger/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
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
  Double_t ti,ti_smeared, tf, tf_smeared, dt;
  Double_t eta, energy;
  Double_t timeResolution;


  const Double_t c_light = 2.99792458E8;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;

    // convert mm in seconds
    ti = candidate->InitialPosition.T()*1.0E-3/c_light;
    tf = candidate->Position.T()*1.0E-3/c_light;

    eta = candidatePosition.Eta();
    energy = candidateMomentum.E();
    timeResolution = fResolutionFormula->Eval(0.0, eta, 0.0, energy);

    dt = timeResolution*gRandom->Gaus(0, 1);
    tf_smeared = tf + dt;
    ti_smeared = ti + dt;
    mother = candidate;
    candidate = static_cast<Candidate*>(candidate->Clone());
    candidate->AddCandidate(mother);

    candidate->Position.SetT(tf_smeared*1.0E3*c_light);
    candidate->ErrorT = timeResolution*1.0E3*c_light;

    // treating charged and neutral differently:
    // for charged we smear the time after propagation, and put a dummy value for time at Vertex
    // since the correct value will be computed after Vertexing4D
    if(candidate->Charge != 0)
      // Dummy Value, correct value will be computed by VertexFinderDA4D
      candidate->InitialPosition.SetT((100+ti)*1.0E3*c_light);
    else
      candidate->InitialPosition.SetT(ti_smeared*1.0E3*c_light);

    fOutputArray->Add(candidate);

  }
}

//------------------------------------------------------------------------------
