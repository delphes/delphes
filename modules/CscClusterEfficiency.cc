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

 /** \class CscClusterEfficiency
  *
  *  This module is specific to the CMS paper searching for neutral LLPs in the CMS endcap muon detectors: https://arxiv.org/abs/2107.04838
  *  It is implemented based on the ClusterEfficiency parameterization function provided in the HEPData entry of the paper: https://www.hepdata.net/record/104408
  *
  *  \author Christina Wang
  *
  */

#include "modules/CscClusterEfficiency.h"

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

CscClusterEfficiency::CscClusterEfficiency() :
  fFormula(0), fItInputArray(0)
{
  fFormula = new DelphesCscClusterFormula;
}

//------------------------------------------------------------------------------

CscClusterEfficiency::~CscClusterEfficiency()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void CscClusterEfficiency::Init()
{
  // read CscClusterEfficiency formula
  fFormula->Compile(GetString("EfficiencyFormula", "1.0"));

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void CscClusterEfficiency::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void CscClusterEfficiency::Process()
{
  Candidate *candidate;
  Double_t Ehad, Eem, decayR, decayZ;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    const TLorentzVector &candidateDecayPosition = candidate->DecayPosition;
    decayZ = abs(candidateDecayPosition.Z());
    decayR = sqrt(pow(candidateDecayPosition.X(),2)+pow(candidateDecayPosition.Y(),2));
    Ehad = candidate->Ehad;
    Eem = candidate->Eem;
    // apply an efficency formula
    if(gRandom->Uniform() > fFormula->Eval(decayR, decayZ, Ehad, Eem)) continue;


    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
