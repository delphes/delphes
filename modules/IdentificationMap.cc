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


/** \class IdentificationMap
 *
 *  Converts particles with some PDG code into another particle,
 *  according to parametrized probability.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/IdentificationMap.h"

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

IdentificationMap::IdentificationMap() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

IdentificationMap::~IdentificationMap()
{
}

//------------------------------------------------------------------------------

void IdentificationMap::Init()
{
  // read efficiency formula


  TMisIDMap::iterator itEfficiencyMap;
  ExRootConfParam param;
  DelphesFormula *formula;
  Int_t i, size, pdg;

  // read efficiency formulas
  param = GetParam("EfficiencyFormula");
  size = param.GetSize();

  fEfficiencyMap.clear();
  for(i = 0; i < size/3; ++i)
  {
    formula = new DelphesFormula;
    formula->Compile(param[i*3 + 2].GetString());
    pdg = param[i*3].GetInt();
    fEfficiencyMap.insert(make_pair(pdg,make_pair(param[i*3 + 1].GetInt(),formula)));

   // cout<<param[i*3].GetInt()<<","<<param[i*3+1].GetInt()<<","<<param[i*3 + 2].GetString()<<endl;

  }

  // set default efficiency formula
  itEfficiencyMap = fEfficiencyMap.find(0);
  if(itEfficiencyMap == fEfficiencyMap.end())
  {
    formula = new DelphesFormula;
    formula->Compile("1.0");

    fEfficiencyMap.insert(make_pair(0,make_pair(0,formula)));
  }

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void IdentificationMap::Finish()
{
  if(fItInputArray) delete fItInputArray;

  TMisIDMap::iterator itEfficiencyMap;
  DelphesFormula *formula;
  for(itEfficiencyMap = fEfficiencyMap.begin(); itEfficiencyMap != fEfficiencyMap.end(); ++itEfficiencyMap)
  {
    formula = (itEfficiencyMap->second).second;
    if(formula) delete formula;
  }

}

//------------------------------------------------------------------------------

void IdentificationMap::Process()
{
  Candidate *candidate;
  Double_t pt, eta, phi;
  TMisIDMap::iterator itEfficiencyMap;
  pair <TMisIDMap::iterator, TMisIDMap::iterator> range;
  DelphesFormula *formula;
  Int_t pdgIn, pdgOut, charge;

  Double_t P, Pi;

 // cout<<"------------ New Event ------------"<<endl;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();
    pt = candidateMomentum.Pt();
    pdgIn = candidate->PID;
    charge = candidate->Charge;

   // cout<<"------------ New Candidate ------------"<<endl;
   // cout<<candidate->PID<<"   "<<pt<<","<<eta<<","<<phi<<endl;

    P = 1.0;

    //first check that PID of this particle is specified in cfg, if not set look for PID=0

    itEfficiencyMap = fEfficiencyMap.find(pdgIn);

    range = fEfficiencyMap.equal_range(pdgIn);
    if(range.first == range.second) range = fEfficiencyMap.equal_range(-pdgIn);
    if(range.first == range.second) range = fEfficiencyMap.equal_range(0);

    //loop over submap for this pid
    for (TMisIDMap::iterator it=range.first; it!=range.second; ++it)
    {

      formula = (it->second).second;
      pdgOut = (it->second).first;

      Pi = formula->Eval(pt, eta);

      // check that sum of probabilities does not exceed 1.
      P = (P - Pi)/P;

      if( P < 0.0 ) continue;
      else
      {

       //randomly assign a PID to particle according to map
       Double_t rndm = gRandom->Uniform();

       if(rndm > P)
       {
         //change PID of particle
         if(pdgOut != 0) candidate->PID = charge*pdgOut;
         fOutputArray->Add(candidate);
         break;
       }
      }

    }

   }

}

//------------------------------------------------------------------------------
