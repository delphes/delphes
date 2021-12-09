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

/** \claCscClusterIdncy
 *
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/CscClusterId.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesCscClusterFormula.h"

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

CscClusterId::CscClusterId() :
  fFormula(0), fItInputArray(0)
{
  fFormula = new DelphesCscClusterFormula;
}

//------------------------------------------------------------------------------

CscClusterId::~CscClusterId()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void CscClusterId::Init()
{
  // read efficiency formula

  fFormula->Compile(GetString("EfficiencyFormula", "1.0"));

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void CscClusterId::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void CscClusterId::Process()
{
  Candidate *candidate;
  Double_t Ehad, Eem, decayR, decayZ, NStationEff, eta, eta_cut;
  Int_t avgStation;
  Double_t signPz, cosTheta;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &candidateDecayPosition = candidate->DecayPosition;
    decayZ = abs(candidateDecayPosition.Z());
    decayR = sqrt(pow(candidateDecayPosition.X(),2)+pow(candidateDecayPosition.Y(),2));
    Ehad = candidate->Ehad;
    Eem = candidate->Eem;

    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signPz * 999.9 : momentum.Eta());

    NStationEff = fFormula->Eval(decayR, decayZ, Ehad); //pt is used as argument in DelphesCscClusterFormula

    // assign average station for the cluster
    if (decayZ < 6320) avgStation = 1;
    else if (decayZ < 7240 && decayR > 2750)avgStation = 1;
    else if (decayZ < 8500) avgStation = 2;
    else if (decayZ < 9700) avgStation = 3;
    else avgStation = 4;

    // if NStation == 1, different eta cut is applied
    if (avgStation == 1) eta_cut = 1.8;
    else if (avgStation == 2) eta_cut = 1.6;
    else if (avgStation == 3) eta_cut = 1.6;
    else if (avgStation == 4) eta_cut = 1.8;
    if(gRandom->Uniform() > NStationEff*(abs(eta)<1.9)+(1.0-NStationEff)*(abs(eta)<eta_cut)) continue;

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
