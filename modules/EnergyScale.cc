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

/** \class EnergyScale
 *
 *  Applies energy scale.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/EnergyScale.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

EnergyScale::EnergyScale() :
  fFormula(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

EnergyScale::~EnergyScale()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void EnergyScale::Init()
{
  // read resolution formula

  fFormula->Compile(GetString("ScaleFormula", "0.0"));

  // import input array(s)
  ImportArray(GetString("InputArray", "FastJetFinder/jets"), fInputArray);
  // create output arrays
  ExportArray(fOutputArray, GetString("OutputArray", "jets"));
}

//------------------------------------------------------------------------------

void EnergyScale::Finish()
{
}

//------------------------------------------------------------------------------

void EnergyScale::Process()
{
  Double_t scale;

  fOutputArray->clear();
  for(const auto &candidate : *fInputArray)
  {
    auto momentum = candidate.Momentum;

    scale = fFormula->Eval(momentum.Pt(), momentum.Eta(), momentum.Phi(), momentum.E());

    if(scale > 0.0) momentum *= scale;

    auto new_candidate = candidate;
    new_candidate.Momentum = momentum;
    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------
