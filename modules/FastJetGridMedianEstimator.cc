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


/** \class FastJetGridMedianEstimator
 *
 *  Computes median energy density per event using a fixed grid.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/FastJetGridMedianEstimator.h"

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
#include <vector>
#include <utility>

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/RectangularGrid.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "fastjet/tools/GridMedianBackgroundEstimator.hh"

#include "fastjet/plugins/SISCone/fastjet/SISConePlugin.hh"
#include "fastjet/plugins/CDFCones/fastjet/CDFMidPointPlugin.hh"
#include "fastjet/plugins/CDFCones/fastjet/CDFJetCluPlugin.hh"

#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"
#include "fastjet/contribs/Nsubjettiness/Njettiness.hh"
#include "fastjet/contribs/Nsubjettiness/NjettinessPlugin.hh"
#include "fastjet/contribs/Nsubjettiness/ExtraRecombiners.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;


//------------------------------------------------------------------------------

FastJetGridMedianEstimator::FastJetGridMedianEstimator() :
  fItInputArray(0)
{

}

//------------------------------------------------------------------------------

FastJetGridMedianEstimator::~FastJetGridMedianEstimator()
{

}

//------------------------------------------------------------------------------

void FastJetGridMedianEstimator::Init()
{
  ExRootConfParam param;
  Long_t i, size;
  Double_t drap, dphi, rapMin, rapMax;

  // read rapidity ranges

  param = GetParam("GridRange");
  size = param.GetSize();

  fEstimators.clear();
  for(i = 0; i < size/4; ++i)
  {
    rapMin = param[i*4].GetDouble();
    rapMax = param[i*4 + 1].GetDouble();
    drap = param[i*4 + 2].GetDouble();
    dphi = param[i*4 + 3].GetDouble();
    fEstimators.push_back(new GridMedianBackgroundEstimator(rapMin, rapMax, drap, dphi));
  }

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "Calorimeter/towers"));
  fItInputArray = fInputArray->MakeIterator();

  fRhoOutputArray = ExportArray(GetString("RhoOutputArray", "rho"));
}

//------------------------------------------------------------------------------

void FastJetGridMedianEstimator::Finish()
{
  vector< GridMedianBackgroundEstimator * >::iterator itEstimators;

  for(itEstimators = fEstimators.begin(); itEstimators != fEstimators.end(); ++itEstimators)
  {
    if(*itEstimators) delete *itEstimators;
  }

  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void FastJetGridMedianEstimator::Process()
{
  Candidate *candidate;
  TLorentzVector momentum;
  Int_t number;
  Double_t rho = 0;
  PseudoJet jet;
  vector< PseudoJet > inputList, outputList;

  vector< GridMedianBackgroundEstimator * >::iterator itEstimators;;

  DelphesFactory *factory = GetFactory();

  inputList.clear();

  // loop over input objects
  fItInputArray->Reset();
  number = 0;
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    momentum = candidate->Momentum;
    jet = PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
    jet.set_user_index(number);
    inputList.push_back(jet);
    ++number;
  }

  // compute rho and store it

  for(itEstimators = fEstimators.begin(); itEstimators != fEstimators.end(); ++itEstimators)
  {
    (*itEstimators)->set_particles(inputList);

    rho = (*itEstimators)->rho();

    candidate = factory->NewCandidate();
    candidate->Momentum.SetPtEtaPhiE(rho, 0.0, 0.0, rho);
    candidate->Edges[0] = (*itEstimators)->rapmin();
    candidate->Edges[1] = (*itEstimators)->rapmax();
    fRhoOutputArray->Add(candidate);
  }
}
