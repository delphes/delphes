/** \class TrackCountingTauTagging
 *
 *  Determines origin of jet,
 *  applies b-tagging efficiency (miss identification rate) formulas
 *  and sets b-tagging flags
 *
 *  $Date: 2013-02-22 01:01:36 +0100 (Fri, 22 Feb 2013) $
 *  $Revision: 926 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFilter.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include "ExRootAnalysis/ExRootClassifier.h"

#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>

using namespace std;

//------------------------------------------------------------------------------

class TrackCountingTauTaggingPartonClassifier: public ExRootClassifier
{
public:
  explicit TrackCountingTauTaggingPartonClassifier(const CandidatesCollection &array) : fParticleInputArray(array) {}

  Int_t GetCategory(TObject *object);

  Double_t fEtaMax, fPTMin;

  const CandidatesCollection &fParticleInputArray;
};

//------------------------------------------------------------------------------

class TrackCountingTauTagging: public DelphesModule
{
public:
  TrackCountingTauTagging() = default;

  void Init() override;
  void Process() override;
  void Finish() override;

private:
  Int_t fBitNumber;

  Double_t fDeltaR;
  Double_t fDeltaRTrack;
  Double_t fTrackPTMin;

  std::map<Int_t, std::unique_ptr<DelphesFormula> > fEfficiencyMap; //!

  std::unique_ptr<TrackCountingTauTaggingPartonClassifier> fClassifier; //!
  std::unique_ptr<DelphesFilter> fFilter;

  CandidatesCollection fParticleInputArray; //!
  CandidatesCollection fTrackInputArray; //!
  CandidatesCollection fPartonInputArray; //!
  CandidatesCollection fJetInputArray; //!
};

//------------------------------------------------------------------------------

Int_t TrackCountingTauTaggingPartonClassifier::GetCategory(TObject *object)
{
  Candidate *tau = static_cast<Candidate *>(object);

  const TLorentzVector &momentum = tau->Momentum;
  Int_t pdgCode, i, j;

  pdgCode = TMath::Abs(tau->PID);
  if(pdgCode != 15) return -1;

  if(momentum.Pt() <= fPTMin || TMath::Abs(momentum.Eta()) > fEtaMax) return -1;

  if(tau->D1 < 0) return -1;

  if(tau->D2 < tau->D1) return -1;

  if(tau->D1 >= static_cast<int>(fParticleInputArray->size()) || tau->D2 >= static_cast<int>(fParticleInputArray->size()))
  {
    throw runtime_error("tau's daughter index is greater than the ParticleInputArray size");
  }

  for(i = tau->D1; i <= tau->D2; ++i)
  {
    Candidate *daughter1 = static_cast<Candidate *>(fParticleInputArray->at(i));
    pdgCode = TMath::Abs(daughter1->PID);
    if(pdgCode == 11 || pdgCode == 13 || pdgCode == 15)
      return -1;
    else if(pdgCode == 24)
    {
      if(daughter1->D1 < 0) return -1;
      for(j = daughter1->D1; j <= daughter1->D2; ++j)
      {
        Candidate *daughter2 = static_cast<Candidate *>(fParticleInputArray->at(j));
        pdgCode = TMath::Abs(daughter2->PID);
        if(pdgCode == 11 || pdgCode == 13) return -1;
      }
    }
  }

  return 0;
}

//------------------------------------------------------------------------------

void TrackCountingTauTagging::Init()
{
  map<Int_t, std::unique_ptr<DelphesFormula> >::iterator itEfficiencyMap;
  ExRootConfParam param;
  Int_t i, size;

  fBitNumber = GetInt("BitNumber", 0);

  fDeltaR = GetDouble("DeltaR", 0.5);
  fDeltaRTrack = GetDouble("DeltaRTrack", 0.2);
  fTrackPTMin = GetDouble("TrackPTMin", 1.0);

  // read efficiency formulas
  param = GetParam("EfficiencyFormula");
  size = param.GetSize();

  fEfficiencyMap.clear();
  for(i = 0; i < size / 2; ++i)
  {
    std::unique_ptr<DelphesFormula> formula = std::make_unique<DelphesFormula>();
    formula->Compile(param[i * 2 + 1].GetString());

    fEfficiencyMap[param[i * 2].GetInt()] = std::move(formula);
  }

  // set default efficiency formula
  itEfficiencyMap = fEfficiencyMap.find(0);
  if(itEfficiencyMap == fEfficiencyMap.end())
  {
    std::unique_ptr<DelphesFormula> formula = std::make_unique<DelphesFormula>();
    formula->Compile("0.0");

    fEfficiencyMap[0] = std::move(formula);
  }

  // import input array(s)
  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "Delphes/allParticles"));
  fPartonInputArray = ImportArray(GetString("PartonInputArray", "Delphes/partons"));
  fTrackInputArray = ImportArray(GetString("TrackInputArray", "TrackMerger/tracks"));
  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));

  fClassifier = std::make_unique<TrackCountingTauTaggingPartonClassifier>(fParticleInputArray);
  fClassifier->fPTMin = GetDouble("TauPTMin", 1.0);
  fClassifier->fEtaMax = GetDouble("TauEtaMax", 2.5);

  fFilter = std::make_unique<DelphesFilter>(fPartonInputArray);
}

//------------------------------------------------------------------------------

void TrackCountingTauTagging::Finish()
{
  fEfficiencyMap.clear();
}

//------------------------------------------------------------------------------

void TrackCountingTauTagging::Process()
{
  TLorentzVector tauMomentum;
  Double_t pt, eta, phi, e;
  map<Int_t, std::unique_ptr<DelphesFormula> >::iterator itEfficiencyMap;
  Int_t charge, i, identifier;

  // select taus
  fFilter->Reset();
  CandidatesCollection tauArray = fFilter->GetSubArray(fClassifier.get(), 0);

  // loop over all input jets
  for(Candidate *const &jet : *fJetInputArray)
  {
    identifier = 0;
    const TLorentzVector &jetMomentum = jet->Momentum;
    charge = 0;
    eta = jetMomentum.Eta();
    phi = jetMomentum.Phi();
    pt = jetMomentum.Pt();
    e = jetMomentum.E();

    // loop over all input tracks
    for(Candidate *const &track : *fTrackInputArray)
    {
      if((track->Momentum).Pt() < fTrackPTMin) continue;
      if(jetMomentum.DeltaR(track->Momentum) <= fDeltaRTrack)
      {
        identifier -= 1;
        charge += track->Charge;
      }
    }

    // loop over all input taus
    bool matchedTau = false;
    for(Candidate *const &tau : *tauArray)
    {
      if(tau->D1 < 0) continue;

      if(tau->D1 >= static_cast<int>(fParticleInputArray->size()) || tau->D2 >= static_cast<int>(fParticleInputArray->size()))
      {
        throw runtime_error("tau's daughter index is greater than the ParticleInputArray size");
      }

      tauMomentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

      for(i = tau->D1; i <= tau->D2; ++i)
      {
        Candidate *daughter = static_cast<Candidate *>(fParticleInputArray->at(i));
        if(TMath::Abs(daughter->PID) == 16) continue;
        tauMomentum += daughter->Momentum;
      }

      if(jetMomentum.DeltaR(tauMomentum) <= fDeltaR)
      {
        matchedTau = true;
      }
    }
    if(matchedTau)
      identifier *= -1;
    // find an efficency formula
    // If the identifier is larger than 2, set it to 2 (multiprong requires at least 2 tracks)
    if(identifier > 2)
      identifier = 2;
    else if(identifier < -2)
      identifier = -2;

    itEfficiencyMap = fEfficiencyMap.find(identifier);
    if(itEfficiencyMap == fEfficiencyMap.end())
    {
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    std::unique_ptr<DelphesFormula> &formula = itEfficiencyMap->second;

    // apply an efficency formula

    // apply an efficency formula
    jet->TauTag |= (gRandom->Uniform() <= formula->Eval(pt, eta, phi, e)) << fBitNumber;

    // set tau charge
    jet->Charge = charge;
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("TrackCountingTauTagging", TrackCountingTauTagging);
