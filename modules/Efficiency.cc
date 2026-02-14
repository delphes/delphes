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

/** \class Efficiency
 *
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Efficiency.h"

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
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

Efficiency::Efficiency() :
  fFormula(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

Efficiency::~Efficiency()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void Efficiency::Init()
{
  // read efficiency formula

  fFormula->Compile(GetString("EfficiencyFormula", "1.0"));

  // import input array(s)
  GetFactory()->EventModel()->Attach(GetString("InputArray", "ParticlePropagator/stableParticles"), fInputArray);

  // switch to compute efficiency based on momentum vector eta, phi
  fUseMomentumVector = GetBool("UseMomentumVector", false);

  // create output array
  GetFactory()->EventModel()->Book(fOutputArray, GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void Efficiency::Finish()
{
}

//------------------------------------------------------------------------------

void Efficiency::Process()
{
  Double_t pt, eta, phi, e;

  for(const auto &candidate : *fInputArray)
  {
    const auto &candidatePosition = candidate.Position;
    const auto &candidateMomentum = candidate.Momentum;
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();

    if(fUseMomentumVector)
    {
      eta = candidateMomentum.Eta();
      phi = candidateMomentum.Phi();
    }

    pt = candidateMomentum.Pt();
    e = candidateMomentum.E();

    // apply an efficency formula
    if(gRandom->Uniform() > fFormula->Eval(pt, eta, phi, e, const_cast<Candidate *>(&candidate))) continue; //TODO: ensure const-qualification

    fOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------
