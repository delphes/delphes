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

/** \class TauTaggingPartonClassifier
 *
 *  Determines jet category
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/TauTaggingPartonClassifier.h"

#include "classes/DelphesClasses.h"

#include <TLorentzVector.h>
#include <TMath.h>

//------------------------------------------------------------------------------

TauTaggingPartonClassifier::TauTaggingPartonClassifier(const TObjArray *array) :
  fParticleInputArray(array) {}

//------------------------------------------------------------------------------

Int_t TauTaggingPartonClassifier::GetCategory(TObject *object)
{
  Candidate *tau = static_cast<Candidate *>(object);
  Candidate *daughter1 = 0;
  Candidate *daughter2 = 0;

  const TLorentzVector &momentum = tau->Momentum;
  Int_t pdgCode, i, j;

  pdgCode = TMath::Abs(tau->PID);
  if(pdgCode != 15) return -1;

  if(momentum.Pt() <= fPTMin || TMath::Abs(momentum.Eta()) > fEtaMax) return -1;

  if(tau->D1 < 0) return -1;

  if(tau->D2 < tau->D1) return -1;

  if(tau->D1 >= fParticleInputArray->GetEntriesFast() || tau->D2 >= fParticleInputArray->GetEntriesFast())
  {
    throw std::runtime_error("tau's daughter index is greater than the ParticleInputArray size");
  }

  for(i = tau->D1; i <= tau->D2; ++i)
  {
    daughter1 = static_cast<Candidate *>(fParticleInputArray->At(i));
    pdgCode = TMath::Abs(daughter1->PID);
    //if(pdgCode == 11 || pdgCode == 13 || pdgCode == 15)
    //  return -1;
    if(pdgCode == 24)
    {
      if(daughter1->D1 < 0) return -1;
      for(j = daughter1->D1; j <= daughter1->D2; ++j)
      {
        daughter2 = static_cast<Candidate *>(fParticleInputArray->At(j));
        pdgCode = TMath::Abs(daughter2->PID);
        if(pdgCode == 11 || pdgCode == 13) return -1;
      }
    }
  }
  return 0;
}
