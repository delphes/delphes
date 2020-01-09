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

/** \class TimeSmearing
 *
 *  Performs transverse momentum resolution smearing.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/TimeSmearing.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

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

TimeSmearing::TimeSmearing() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

TimeSmearing::~TimeSmearing()
{  
}

//------------------------------------------------------------------------------

void TimeSmearing::Init()
{
  // read resolution formula

  fTimeResolution = GetDouble("TimeResolution", 1.);

  // import input array
  fEtaMax = GetDouble("EtaMax", 6.);
  fInputArray = ImportArray(GetString("InputArray", "MuonMomentumSmearing/muons"));
  fItInputArray = fInputArray->MakeIterator();
  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "muons"));
}

//------------------------------------------------------------------------------

void TimeSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void TimeSmearing::Process()
{
  Candidate *candidate, *mother;
  Double_t ti, tf_smeared, tf;
  Double_t pt, eta, phi, e, p, l;
  Double_t sigma_t, beta_particle;


  const Double_t c_light = 2.99792458E8;

  cout << " STARTINNNGG ---------->" << endl;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    ti = candidate->InitialPosition.T()*1.0E-3/c_light;
    tf = candidate->Position.T()*1.0E-3/c_light;
    
    // dummy, only need to properly call TFormula
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();
    pt = candidateMomentum.Pt();
    e = candidateMomentum.E();
    p = candidateMomentum.P();
    beta_particle = p/e;
    l = candidate->L;
    // apply smearing formula

    if(candidate->Charge != 0)
    { 
      tf_smeared = tf + fTimeResolution*gRandom->Gaus(0, 1);
      //mother = candidate;
      //candidate = static_cast<Candidate*>(candidate->Clone());  // I am not sure that we need these lines !!!
      //candidate->AddCandidate(mother);
      candidate->InitialPosition.SetT((100+ti)*1.0E3*c_light);
      candidate->Position.SetT(tf_smeared*1.0E3*c_light);
      candidate->ErrorT = fTimeResolution*1.0E3*c_light;
      fOutputArray->Add(candidate);
    }
    else
    {
    	sigma_t = sqrt(pow(20,2) + pow(150,2)/e);
    	ti = sigma_t - l*1.0E3/(c_light*beta_particle);   
    	candidate->InitialPosition.SetT(ti);
    	candidate->ErrorT = sigma_t*1.0E3*c_light;		// Do we need to sum with 100 like in upside ?
    	fOutputArray->Add(candidate);
    }
   
  }
}

//------------------------------------------------------------------------------
