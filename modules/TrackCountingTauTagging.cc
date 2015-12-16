
/** \class TrackCountingTauTagging
 *
 *  Determines origin of jet,
 *  applies b-tagging efficiency (miss identification rate) formulas
 *  and sets b-tagging flags
 *
 *  $Date: 2013-07-12 00:22:27 +0200 (Fri, 12 Jul 2013) $
 *  $Revision: 1217 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/TrackCountingTauTagging.h"

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

class TrackCountingTauTaggingPartonClassifier : public ExRootClassifier
{
public:

  TrackCountingTauTaggingPartonClassifier(const TObjArray *array);

  Int_t GetCategory(TObject *object);

  Double_t fEtaMax, fPTMin;

  const TObjArray *fParticleInputArray;
};

//------------------------------------------------------------------------------
TrackCountingTauTaggingPartonClassifier::TrackCountingTauTaggingPartonClassifier(const TObjArray *array) :
  fParticleInputArray(array)
{
}

//------------------------------------------------------------------------------

Int_t TrackCountingTauTaggingPartonClassifier::GetCategory(TObject *object)
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

  if(tau->D1 >= fParticleInputArray->GetEntriesFast() ||
     tau->D2 >= fParticleInputArray->GetEntriesFast())
  {
    throw runtime_error("tau's daughter index is greater than the ParticleInputArray size");
  }

  for(i = tau->D1; i <= tau->D2; ++i)
  {
    daughter1 = static_cast<Candidate *>(fParticleInputArray->At(i));
    pdgCode = TMath::Abs(daughter1->PID);
    if(pdgCode == 11 || pdgCode == 13 || pdgCode == 15) return -1;
    else if(pdgCode == 24)
    {
     if(daughter1->D1 < 0) return -1;
     for(j = daughter1->D1; j <= daughter1->D2; ++j)
     {
       daughter2 = static_cast<Candidate*>(fParticleInputArray->At(j));
       pdgCode = TMath::Abs(daughter2->PID);
       if(pdgCode == 11 || pdgCode == 13) return -1;
     }
	
    }
  }

  return 0;
}

//------------------------------------------------------------------------------

TrackCountingTauTagging::TrackCountingTauTagging() :
  fClassifier(0), fFilter(0),
  fItPartonInputArray(0), fItTrackInputArray(0), fItJetInputArray(0)
{
}

//------------------------------------------------------------------------------

TrackCountingTauTagging::~TrackCountingTauTagging()
{
}

//------------------------------------------------------------------------------

void TrackCountingTauTagging::Init()
{
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  ExRootConfParam param;
  DelphesFormula *formula;
  Int_t i, size;

  fBitNumber = GetInt("BitNumber", 0);

  fDeltaR = GetDouble("DeltaR", 0.5);
  fDeltaRTrack = GetDouble("DeltaRTrack", 0.2);
  fTrackPTMin = GetDouble("TrackPTMin", 1.0);
  
  // read efficiency formulas
  param = GetParam("EfficiencyFormula");
  size = param.GetSize();

  fEfficiencyMap.clear();
  for(i = 0; i < size/2; ++i)
  {
    formula = new DelphesFormula;
    formula->Compile(param[i*2 + 1].GetString());

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

  // import input array(s)

  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "Delphes/allParticles"));

  fClassifier = new TrackCountingTauTaggingPartonClassifier(fParticleInputArray);
  fClassifier->fPTMin = GetDouble("TauPTMin", 1.0);
  fClassifier->fEtaMax = GetDouble("TauEtaMax", 2.5);

  fPartonInputArray = ImportArray(GetString("PartonInputArray", "Delphes/partons"));
  fItPartonInputArray = fPartonInputArray->MakeIterator();

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "TrackMerger/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();
  
  fFilter = new ExRootFilter(fPartonInputArray);

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();
}

//------------------------------------------------------------------------------

void TrackCountingTauTagging::Finish()
{
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  DelphesFormula *formula;

  if(fFilter) delete fFilter;
  if(fClassifier) delete fClassifier;
  if(fItJetInputArray) delete fItJetInputArray;
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItPartonInputArray) delete fItPartonInputArray;

  for(itEfficiencyMap = fEfficiencyMap.begin(); itEfficiencyMap != fEfficiencyMap.end(); ++itEfficiencyMap)
  {
    formula = itEfficiencyMap->second;
    if(formula) delete formula;
  }
}

//------------------------------------------------------------------------------

void TrackCountingTauTagging::Process()
{
  Candidate *jet, *tau, *track, *daughter;
  TLorentzVector tauMomentum;
  Double_t pt, eta, phi, e;
  TObjArray *tauArray;
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  DelphesFormula *formula;
  Int_t pdgCode, charge, i, identifier;

  // select taus
  fFilter->Reset();
  tauArray = fFilter->GetSubArray(fClassifier, 0);

  if(tauArray == 0) return;

  TIter itTauArray(tauArray);

  // loop over all input jets
  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate *>(fItJetInputArray->Next())))
  {
    identifier = 0;
    const TLorentzVector &jetMomentum = jet->Momentum;
    pdgCode = 0;
    charge = 0;
    eta = jetMomentum.Eta();
    phi = jetMomentum.Phi();
    pt = jetMomentum.Pt();
    e = jetMomentum.E();


// loop over all input tracks
    fItTrackInputArray->Reset();
    while((track = static_cast<Candidate *>(fItTrackInputArray->Next())))
    {
        if((track->Momentum).Pt() < fTrackPTMin) continue;
        if(jetMomentum.DeltaR(track->Momentum) <= fDeltaRTrack) {
            identifier -= 1;
            charge += track->Charge;
        }
    }
    
    // loop over all input taus
    itTauArray.Reset();
    bool matchedTau = false;
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

      if(jetMomentum.DeltaR(tauMomentum) <= fDeltaR)
      {        
        matchedTau = true;
        pdgCode = 15;
      }
    }
    if(matchedTau)
	identifier *= -1;
    // find an efficency formula
    // If the identifier is larger than 2, set it to 2 (multiprong requires at least 2 tracks)
    if (identifier > 2)
	identifier = 2;
    else if (identifier < -2)
	identifier = -2;

    
    itEfficiencyMap = fEfficiencyMap.find(identifier);
    if(itEfficiencyMap == fEfficiencyMap.end())
    {
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    formula = itEfficiencyMap->second;

    // apply an efficency formula

    // apply an efficency formula
    jet->TauTag |= (gRandom->Uniform() <= formula->Eval(pt, eta, phi, e)) << fBitNumber;
   
   
    // set tau charge
    jet->Charge = charge;
  }
}

//------------------------------------------------------------------------------
