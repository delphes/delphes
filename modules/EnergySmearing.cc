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


/** \class EnergySmearing
 *
 *  Performs energy resolution smearing.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/EnergySmearing.h"

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

EnergySmearing::EnergySmearing() :
  fFormula(0), fItInputArray(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

EnergySmearing::~EnergySmearing()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void EnergySmearing::Init()
{
  // read resolution formula

  fFormula->Compile(GetString("ResolutionFormula", "0.0"));

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void EnergySmearing::Finish()
{  
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void EnergySmearing::Process()
{
  Candidate *candidate, *mother;
  Double_t pt, energy, eta, phi;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    
    pt = candidatePosition.Pt();
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();
    energy = candidateMomentum.E();
 
    // apply smearing formula
    energy = gRandom->Gaus(energy, fFormula->Eval(pt, eta, phi, energy));
     
    if(energy <= 0.0) continue;
 
    mother = candidate;
    candidate = static_cast<Candidate*>(candidate->Clone());
    eta = candidateMomentum.Eta();
    phi = candidateMomentum.Phi();
    candidate->Momentum.SetPtEtaPhiE(energy/TMath::CosH(eta), eta, phi, energy);
    candidate->TrackResolution = fFormula->Eval(pt, eta, phi, energy)/candidateMomentum.E();
    candidate->AddCandidate(mother);
 
    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
