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

/** \class MetaStableCharged
 *
 *  Selects candidates from the InputArray that are consistent with a 
 *  meta-stable charged particle (status != 1) travelling beyond 
 *  a certain cylinder (assumed to collect min Nhits for tracking) 
 *
 *  \author S.Vashishtha, M.Goblirsch-Kolb - DESY
 *
 */

#include "modules/MetaStableCharged.h"

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

MetaStableCharged::MetaStableCharged() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

MetaStableCharged::~MetaStableCharged()
{
}

//------------------------------------------------------------------------------

void MetaStableCharged::Init()
{
  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "metastableParticles"));
  fMinTravel = GetDouble("MinTravel", 1.0);
}

//------------------------------------------------------------------------------

void MetaStableCharged::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void MetaStableCharged::Process()
{
  Candidate *candidate;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    // discard neutrals 
    if (std::abs(candidate->Charge) < 1.e-9) continue;
    // discard stable - already accounted for
    if (candidate->Status == 1) continue;
    // cut on travel distance - using stored decay Length L 
    if ( candidate->L >= fMinTravel*1e3 ){    
      fOutputArray->Add(candidate);
    }
  } 

}
//------------------------------------------------------------------------------
