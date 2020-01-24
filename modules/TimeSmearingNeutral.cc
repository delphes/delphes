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

/** \class TimeSmearingNeutral
 *
 *  Performs transverse momentum resolution smearing.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/TimeSmearingNeutral.h"

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

TimeSmearingNeutral::TimeSmearingNeutral() :
  fFormula(0), fItInputArray(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

TimeSmearingNeutral::~TimeSmearingNeutral()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void TimeSmearingNeutral::Init()
{
  // read resolution formula

  //fFormula->Compile(GetString("TimeResolution", "1.0"));
  fTimeResolution = GetDouble("TimeResolution", 1.);
  // import input array
  fEtaMax = GetDouble("EtaMax", 6.);
  fInputArray = ImportArray(GetString("InputArray", "MuonMomentumSmearing/muons"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "muons"));
}

//------------------------------------------------------------------------------

void TimeSmearingNeutral::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void TimeSmearingNeutral::Process()
{
  Candidate *candidate, *mother;
  Double_t ti, tf_smeared, tf, timeResolution;
  Double_t pt, eta, phi, e, p, d0, dz, ctgTheta;
  Double_t beta_particle;
  Double_t l;
  Double_t sigma_t;

  const Double_t c_light = 2.99792458E8;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    ti = candidate->InitialPosition.T()*1.0E-3/c_light;
    tf = candidate->Position.T()*1.0E-3/c_light;
    e = candidate->Momentum.E();
    p = candidate->Momentum.P();
    beta_particle = p/e;
    l = candidate->L;

    if(candidate->Charge == 0)
    {
      sigma_t = sqrt(pow(20,2) + pow(150,2)/e)*1.0E-12;
    }else
      sigma_t = fTimeResolution*gRandom->Gaus(0, 1);;

    ti = sigma_t - l*1.0E3/(c_light*beta_particle);


    candidate->InitialPosition.SetT(ti);
    candidate->ErrorT = sigma_t*1.0E3*c_light;
    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
