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

/** \class DecayFilter
 *
 *  This module randomly generates decays along the particle trajectory length
 *  according to actual particle decay length, taking into account for the boost
 *  and using ROOT TDatabasePDG as a source for the particle lifetime.
 *
 *  This module is to to be used after a PropagateParticle step or a similar module
 *  that calculates and store a trajectory length.
 *
 *  Particles that decay are not added to the OutputArray.
 *
 *  \author R. Preghenella - INFN, Bologna
 *
 */

#include "modules/DecayFilter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TParticlePDG.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

void DecayFilter::Init()
{
  // import input array(s)
  GetFactory()->EventModel()->Attach(GetString("InputArray", "FastJetFinder/jets"), fInputArray);
  // create output arrays
  GetFactory()->EventModel()->Book(fOutputArray, GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void DecayFilter::Finish()
{
}

//------------------------------------------------------------------------------

void DecayFilter::Process()
{
  TDatabasePDG *pdgdb = TDatabasePDG::Instance();
  const Double_t c = TMath::C(); // [m/s]
  Double_t m, t, p, bgct, L, l;

  // loop over all input candidates
  for(const auto &candidate : *fInputArray)
  {
    // get particle information from PDG
    TParticlePDG *pdg = pdgdb->GetParticle(candidate.PID);
    if(!pdg)
    { // don't know this particle
      fOutputArray->emplace_back(candidate);
      continue;
    }
    m = pdg->Mass();
    t = pdg->Lifetime(); // [s]
    if(t == 0.)
    { // does not decay
      fOutputArray->emplace_back(candidate);
      continue;
    }

    // compute boosted decay length (beta gamma c tau)
    p = candidate.P;
    bgct = p / m * c * t; // [m]

    // get full trajectory length and generate random decay length
    L = candidate.L * 1.0E-3; // [m]
    l = gRandom->Exp(bgct);

    // if random decay happens before end of trajectory, reject track
    if(l < L) continue;

    // else particle did not decay within the trajectory
    fOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------
