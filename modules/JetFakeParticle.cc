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


/** \class JetFakeParticle
 *
 *  Converts jet into particle with some PID,
 *  according to parametrized probability.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/JetFakeParticle.h"

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

JetFakeParticle::JetFakeParticle() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

JetFakeParticle::~JetFakeParticle()
{
}

//------------------------------------------------------------------------------

void JetFakeParticle::Init()
{
  TFakeMap::iterator itEfficiencyMap;
  ExRootConfParam param;
  DelphesFormula *formula;
  Int_t i, size, pdgCode;

  // read efficiency formulas
  param = GetParam("EfficiencyFormula");
  size = param.GetSize();

  fEfficiencyMap.clear();

  for(i = 0; i < size/2; ++i)
  {
    formula = new DelphesFormula;
    formula->Compile(param[i*2 + 1].GetString());
    pdgCode = param[i*2].GetInt();

    if(TMath::Abs(pdgCode) != 11 && TMath::Abs(pdgCode) != 13 && TMath::Abs(pdgCode) != 22)
    {
      throw runtime_error("Jets can only fake into electrons, muons or photons. Other particles are not authorized.");
    }

    fEfficiencyMap[param[i*2].GetInt()] = formula;
  }

  // set default efficiency formula
  itEfficiencyMap = fEfficiencyMap.find(0);
  if(itEfficiencyMap == fEfficiencyMap.end())
  {
    formula = new DelphesFormula;
    formula->Compile("0.0");

    fEfficiencyMap[0] = formula;
  }

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "FastJetFinder/jets"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fElectronOutputArray = ExportArray(GetString("ElectronOutputArray", "fakeElectrons"));
  fMuonOutputArray = ExportArray(GetString("MuonOutputArray", "fakeMuons"));
  fPhotonOutputArray = ExportArray(GetString("PhotonOutputArray", "fakePhotons"));
  fJetOutputArray = ExportArray(GetString("JetOutputArray", "jets"));

}

//------------------------------------------------------------------------------

void JetFakeParticle::Finish()
{
  if(fItInputArray) delete fItInputArray;

  TFakeMap::iterator itEfficiencyMap;
  DelphesFormula *formula;
  for(itEfficiencyMap = fEfficiencyMap.begin(); itEfficiencyMap != fEfficiencyMap.end(); ++itEfficiencyMap)
  {
    formula = itEfficiencyMap->second;
    if(formula) delete formula;
  }
}

//------------------------------------------------------------------------------

void JetFakeParticle::Process()
{
  Candidate *candidate, *fake = 0;
  Double_t pt, eta, phi, e;
  TFakeMap::iterator itEfficiencyMap;
  DelphesFormula *formula;
  Int_t pdgCodeOut;

  Double_t p, r, rs, total;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();
    pt = candidateMomentum.Pt();
    e = candidateMomentum.E();

    r = gRandom->Uniform();
    total = 0.0;
    fake = 0;

    // loop over map for this jet
    for(itEfficiencyMap = fEfficiencyMap.begin(); itEfficiencyMap != fEfficiencyMap.end(); ++itEfficiencyMap)
    {
      formula = itEfficiencyMap->second;
      pdgCodeOut = itEfficiencyMap->first;

      p = formula->Eval(pt, eta, phi, e);

      if(total <= r && r < total + p)
      {
        fake = static_cast<Candidate*>(candidate->Clone());

        // convert jet

        if(TMath::Abs(pdgCodeOut) == 11 || TMath::Abs(pdgCodeOut) == 13)
        {
          if(candidate->Charge != 0)
          {
            fake->Charge = candidate->Charge/TMath::Abs(candidate->Charge);
          }
          else
          {
            rs = gRandom->Uniform();
            fake->Charge = (rs < 0.5) ? -1 : 1;
            
          }
        }

        if(TMath::Abs(pdgCodeOut) == 22) fake->PID = 22;

        if(TMath::Abs(pdgCodeOut) == 11) fElectronOutputArray->Add(fake);
        if(TMath::Abs(pdgCodeOut) == 13) fMuonOutputArray->Add(fake);
        if(TMath::Abs(pdgCodeOut) == 22) fPhotonOutputArray->Add(fake);

        break;
      }

      total += p;
    }

    if(!fake) fJetOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
