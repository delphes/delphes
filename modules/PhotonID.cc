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

/** \class PhotonID
 *
 *  Applies complex photon Id. Reconstructed photon candidtes are first separated into matched and non-matched to gen particles. 
 *  Non-matched pass the "fake" efficiency. Matched photons get further splitted into isolated and non-isolated (user can choose criterion for isolation)
 *  Isolated photons pass the "prompt" efficiency while the non-isolated pass the "non-prompt" efficiency
 *
 *  \author M. Selvaggi CERN
 *
 */

#include "modules/PhotonID.h"

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

PhotonID::PhotonID() :
  fPromptFormula(0), fNonPromptFormula(0), fFakeFormula(0), fItInputPhotonArray(0), fItInputGenArray(0)
{
  fPromptFormula = new DelphesFormula;
  fNonPromptFormula = new DelphesFormula;
  fFakeFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

PhotonID::~PhotonID()
{
  if(fPromptFormula) delete fPromptFormula;
  if(fNonPromptFormula) delete fNonPromptFormula;
  if(fFakeFormula) delete fFakeFormula;
}

//------------------------------------------------------------------------------

void PhotonID::Init()
{

  // read PhotonID formulae
  fPromptFormula->Compile(GetString("PromptFormula", "1.0"));
  fNonPromptFormula->Compile(GetString("NonPromptFormula", "1.0"));
  fFakeFormula->Compile(GetString("FakeFormula", "1.0"));

  // import input arrays
  fInputPhotonArray = ImportArray(GetString("InputPhotonArray", "PhotonIsolation/photons"));
  fItInputPhotonArray = fInputPhotonArray->MakeIterator();

  // use filtered collection for speed
  fInputGenArray = ImportArray(GetString("InputGenArray", "GenParticleFilter/filteredParticles"));
  fItInputGenArray = fInputGenArray->MakeIterator();

  // min pt to be considered, make sure this threshold is higher than threshold in particle filter
  fPTMin = GetDouble("PTMin", 10.0);

  // to be tuned, since FS and delphes have different isolation profiles
  fRelIsoMax = GetDouble("fRelIsoMax", 0.3);

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "photons"));
}

//------------------------------------------------------------------------------

void PhotonID::Finish()
{
  if(fItInputPhotonArray) delete fItInputPhotonArray;
  if(fItInputGenArray) delete fItInputGenArray;
}

//------------------------------------------------------------------------------

void PhotonID::Process()
{
  Candidate *candidate, *mother;
  Double_t pt, eta, phi, e;
  Double_t relIso;
  Bool_t isolated;

  //cout<< "----  new event ---------"<<endl;

  fItInputPhotonArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputPhotonArray->Next())))
  {

    mother = candidate;
    candidate = static_cast<Candidate *>(candidate->Clone());
    candidate->AddCandidate(mother);

    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();
    pt = candidateMomentum.Pt();
    e = candidateMomentum.E();

    if(pt < fPTMin) continue;

    //cout<< "              ---- photon -----: "<<pt<<","<<eta<<","<<phi<<endl;

    // find out if photon matches does not match photon in gen collection and apply fae efficiency
    if(isFake(candidate))
    {
      //cout<<"                    Fake!"<<endl;

      if(gRandom->Uniform() > fFakeFormula->Eval(pt, eta, phi, e)) continue;
      //cout<<"                    passed"<<endl;
      candidate->Status = 3;
      fOutputArray->Add(candidate);
    }

    // if matches photon in gen collection
    else
    {
      relIso = candidate->IsolationVar;
      isolated = (relIso < 0.3);
      //cout<<"                    Prompt!:   "<<relIso<<endl;

      // if isolated apply prompt formula
      if(isolated)
      {
        //cout<<"                       isolated!:   "<<relIso<<endl;
        if(gRandom->Uniform() > fPromptFormula->Eval(pt, eta, phi, e)) continue;
        //cout<<"                       passed"<<endl;
        candidate->Status = 1;
        fOutputArray->Add(candidate);
      }

      // if non-isolated apply non-prompt formula
      else
      {
        //cout<<"                       non-isolated!:   "<<relIso<<endl;
        if(gRandom->Uniform() > fNonPromptFormula->Eval(pt, eta, phi, e)) continue;
        //cout<<"                       passed"<<endl;
        candidate->Status = 2;
        fOutputArray->Add(candidate);
      }
    }
  }
}

//------------------------------------------------------------------------------

Bool_t PhotonID::isFake(const Candidate *obj)
{

  const TLorentzVector &mom_rec = obj->Momentum;

  Bool_t matches = false;
  fItInputGenArray->Reset();
  Candidate *gen;

  while((gen = static_cast<Candidate *>(fItInputGenArray->Next())))
  {
    const TLorentzVector &mom_gen = gen->Momentum;
    Int_t status = gen->Status;
    Int_t pdgCode = TMath::Abs(gen->PID);
    Float_t dPtOverPt = TMath::Abs((mom_gen.Pt() - mom_rec.Pt()) / mom_rec.Pt());
    Float_t deltaR = mom_gen.DeltaR(mom_rec);

    if(status != 1) continue;
    if(pdgCode != 22) continue;
    if(dPtOverPt > 0.5) continue;
    if(deltaR > 0.1) continue;

    matches = true;
    break;
  }

  return !matches;
}
