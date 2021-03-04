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

/** \class TimeSmearingTail
 *
 *  Performs transverse momentum resolution smearing.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 *  With addition of non-gaussian tails (R+Preghenella, preghenella@bo.infn.it)
 *  non-gaussian tails can be added on both the left/right side of the peak
 *  according to the "TailLeft" and "TailRight" parameters
 *  following a q-gaussian PDF, where Tail is the value of q (=1 for gaussian)
 *
 */

#include "modules/TimeSmearingTail.h"

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
#include "TF1.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

TimeSmearingTail::TimeSmearingTail() :
  fSignal(0),
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

TimeSmearingTail::~TimeSmearingTail()
{
  if (fSignal)
    delete fSignal;
}

//------------------------------------------------------------------------------

void TimeSmearingTail::Init()
{
  // read resolution formula

  fTimeResolution = GetDouble("TimeResolution", 1.0E-10);
  fTailLeft = GetDouble("TailLeft", 1.0);
  fTailRight = GetDouble("TailRight", 1.0);
  // import input array

  fInputArray = ImportArray(GetString("InputArray", "MuonMomentumSmearing/muons"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "muons"));

  // create signal function

  fSignal = new TF1("fSignal", TimeSmearingTail::qTOF, -10.*fTimeResolution, 20.*fTimeResolution, 3);
  fSignal->SetParameter(0, fTimeResolution);
  fSignal->SetParameter(1, fTailRight);
  fSignal->SetParameter(2, fTailLeft);
}

//------------------------------------------------------------------------------

void TimeSmearingTail::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void TimeSmearingTail::Process()
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

    tf_smeared = tf + fSignal->GetRandom();
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
