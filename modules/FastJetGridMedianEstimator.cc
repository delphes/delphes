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
#include "fastjet/contribs/Nsubjettiness/WinnerTakeAllRecombiner.hh"

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
  // read eta ranges

  ExRootConfParam param = GetParam("GridRange");
  Long_t i, size;
  
  fGrid.clear();
  size = param.GetSize();
  for(i = 0; i < size/4; ++i)
  {
    fGrid[make_pair(param[i*4].GetDouble(), param[i*4 + 1].GetDouble())] = make_pair(param[i*4 + 2].GetDouble(), param[i*4 + 3].GetDouble());
    //cout<<param[i*4].GetDouble()<<","<<  param[i*4 + 1].GetDouble()<<","<< param[i*4 + 2].GetDouble()<<","<< param[i*4 + 3].GetDouble()<<endl;

  }
 
  
  //cout<<fGrid[make_pair(0.0,2.5)].first<<","<<fGrid[make_pair(0.0,2.5)].second<<endl;
 
 // import input array

  fInputArray = ImportArray(GetString("InputArray", "Calorimeter/towers"));
  fItInputArray = fInputArray->MakeIterator();

  fRhoOutputArray = ExportArray(GetString("RhoOutputArray", "rho"));
  
}

//------------------------------------------------------------------------------

void FastJetGridMedianEstimator::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void FastJetGridMedianEstimator::Process()
{
  Candidate *candidate;
  TLorentzVector momentum;
  Double_t deta, dphi, detaMin, detaMax;
  Int_t number;
  Double_t rho = 0;
  PseudoJet jet;
  vector<PseudoJet> inputList, outputList;
  
  std::map< std::pair< Double_t , Double_t > , std::pair< Double_t , Double_t > >::iterator itGrid;
  
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
  
  // cout<<"caio"<<endl;
  for(itGrid = fGrid.begin(); itGrid != fGrid.end(); ++itGrid)
  {
   //Selector select_rapidity = SelectorAbsRapRange(itEtaRangeMap->first, itEtaRangeMap->second);
   // JetMedianBackgroundEstimator estimator(select_rapidity, *fDefinition, *fAreaDefinition);
      
    //cout<<itGrid->first.first<<endl;
   
    detaMin = (itGrid->first).first;
    detaMax = (itGrid->first).second;
    deta    = (itGrid->second).first;
    dphi    = (itGrid->second).second;
    
    //cout<<detaMin<<","<<detaMax<<","<<deta<<","<<dphi<<endl;
    
    
    RectangularGrid grid(detaMin, detaMax, deta, dphi);
    //cout<<grid.is_initialised()<<endl;
    //cout<<grid.rapmin()<<","<<grid.rapmax()<<","<<grid.drap()<<","<<grid.dphi()<<endl;
      
    GridMedianBackgroundEstimator estimator(grid);
      
    estimator.set_particles(inputList);
      
    //cout<<estimator.description()<<endl;
      
    rho = estimator.rho();
    //cout<<rho<<endl;      


    candidate = factory->NewCandidate();
    candidate->Momentum.SetPtEtaPhiE(rho, 0.0, 0.0, rho);
    candidate->Edges[0] = detaMin;
    candidate->Edges[1] = detaMax;
    fRhoOutputArray->Add(candidate);
  }
 
}
