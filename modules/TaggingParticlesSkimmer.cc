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

//------------------------------------------------------------------------------

/** \class TaggingParticlesSkimmer
 *
 *  Filters particle collection by only keeping gen particles useful for b and tau tagging.
 *  tau particles are replaced by their "visible decay".
 *
 *  \author M. Selvaggi
 *
 */

#include "modules/TaggingParticlesSkimmer.h"
#include "modules/TauTagging.h"

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

TaggingParticlesSkimmer::TaggingParticlesSkimmer() :
  fItPartonInputArray(0), fFilter(0), fClassifier(0)
{
}

//------------------------------------------------------------------------------

TaggingParticlesSkimmer::~TaggingParticlesSkimmer()
{
}

//------------------------------------------------------------------------------

void TaggingParticlesSkimmer::Init()
{
 
  fPTMin = GetDouble("PTMin", 15.0);
  fEtaMax = GetDouble("EtaMax", 2.5);
   
  // import input array
  fPartonInputArray = ImportArray(GetString("PartonInputArray", "Delphes/partons"));
  fItPartonInputArray = fPartonInputArray->MakeIterator();
  
  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "Delphes/allParticles"));

 
  fClassifier = new TauTaggingPartonClassifier(fParticleInputArray);
  fClassifier->fPTMin = GetDouble("PTMin", 15.0);
  fClassifier->fEtaMax = GetDouble("EtaMax", 2.5);

 
  fFilter = new ExRootFilter(fPartonInputArray);

  // output array
  fOutputArray = ExportArray(GetString("OutputArray", "taggingParticles"));
}

//------------------------------------------------------------------------------

void TaggingParticlesSkimmer::Finish()
{ 
  if(fItPartonInputArray) delete fItPartonInputArray;
  if(fFilter) delete fFilter;
  if(fClassifier) delete fClassifier;
}

//------------------------------------------------------------------------------

void TaggingParticlesSkimmer::Process()
{
  Candidate *candidate, *tau, *daughter;
  TLorentzVector tauMomentum;
  Double_t pt, eta;
  TObjArray *tauArray;
  Int_t pdgCode, i;

  // first select hadronic taus and replace them by visible part
  fFilter->Reset();
  tauArray = fFilter->GetSubArray(fClassifier, 0);

  if(tauArray == 0) return;

  TIter itTauArray(tauArray);

  // loop over all input taus
  itTauArray.Reset();
  while((tau = static_cast<Candidate *>(itTauArray.Next())))
  {
    if(tau->D1 < 0) continue;

    if(tau->D1 >= fParticleInputArray->GetEntriesFast() ||
       tau->D2 >= fParticleInputArray->GetEntriesFast())
    {
      throw runtime_error("tau's daughter index is greater than the ParticleInputArray size");
    }

    tauMomentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    
    for(i = tau->D1; i <= tau->D2; ++i)
    {
      daughter = static_cast<Candidate *>(fParticleInputArray->At(i));
      if(TMath::Abs(daughter->PID) == 16) continue;
      tauMomentum += daughter->Momentum;
    }
  
   candidate = static_cast<Candidate*>(tau->Clone());
   candidate->Momentum = tauMomentum;

   
   fOutputArray->Add(candidate);

  }

  // then add all other partons (except tau's to avoid double counting)
  
  fItPartonInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItPartonInputArray->Next())))
  {
    pdgCode = TMath::Abs(candidate->PID);
    if(pdgCode == 15) continue;
   
    pt = candidate->Momentum.Pt();
    if(pt < fPTMin) continue;
   
    eta = TMath::Abs(candidate->Momentum.Eta());
    if(eta > fEtaMax) continue;
        
    fOutputArray->Add(candidate);
  }


}

