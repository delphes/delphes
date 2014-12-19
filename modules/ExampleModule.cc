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


/** \class ExampleModule
 *
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/ExampleModule.h"

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

ExampleModule::ExampleModule() :
  fFormula(0), fItInputArray(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

ExampleModule::~ExampleModule()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void ExampleModule::Init()
{
  // read parameters

  fIntParam = GetInt("IntParam", 10);

  fDoubleParam = GetDouble("DoubleParam", 1.0);

  fFormula->Compile(GetString("EfficiencyFormula", "0.4"));

  ExRootConfParam param = GetParam("ArrayParam");
  Long_t i, size;
  fArrayParam.clear();

  size = param.GetSize();
  for(i = 0; i < size; ++i)
  {
    fArrayParam.push_back(param[i].GetDouble());
  }
  
  // import input array(s)

  fInputArray = ImportArray(GetString("InputArray", "FastJetFinder/jets"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array(s)

  fOutputArray = ExportArray(GetString("OutputArray", "jets"));

}

//------------------------------------------------------------------------------

void ExampleModule::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void ExampleModule::Process()
{
  Candidate *candidate;
  TLorentzVector candidatePosition, candidateMomentum;

  // loop over all input candidates
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    candidatePosition = candidate->Position;
    candidateMomentum = candidate->Momentum;

    // apply an efficency formula
    if(gRandom->Uniform() <= fFormula->Eval(candidateMomentum.Pt(), candidatePosition.Eta()))
    {
      fOutputArray->Add(candidate);
    }
  }
}

//------------------------------------------------------------------------------
