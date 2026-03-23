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
#include <TRandom3.h>

using namespace std;

//------------------------------------------------------------------------------

class TrackCountingTauTaggingPartonClassifier: public ExRootClassifier
{
public:
  explicit TrackCountingTauTaggingPartonClassifier(const CandidatesCollection &array) : fParticleInputArray(array) {}

  int GetCategory(TObject *object);

  double fEtaMax, fPTMin;

  const CandidatesCollection &fParticleInputArray;
};

//------------------------------------------------------------------------------

class TrackCountingTauTagging: public DelphesModule
{
public:
  explicit TrackCountingTauTagging(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fBitNumber(Steer<int>("BitNumber", 0)),
    //
    fDeltaR(Steer<double>("DeltaR", 0.5)),
    fDeltaRTrack(Steer<double>("DeltaRTrack", 0.2)),
    fTrackPTMin(Steer<double>("TrackPTMin", 1.0))
  {
    for(const std::pair<int, std::string> &efficiencyFormula :
      Steer<std::vector<std::pair<int, std::string> > >("EfficiencyFormula"))
    {
      std::unique_ptr<DelphesFormula> formula = std::make_unique<DelphesFormula>();
      formula->Compile(efficiencyFormula.second);
      fEfficiencyMap[efficiencyFormula.first] = std::move(formula);
    }
    // set default efficiency formula
    if(fEfficiencyMap.count(0) == 0)
    {
      std::unique_ptr<DelphesFormula> formula = std::make_unique<DelphesFormula>();
      formula->Compile("0.0");
      fEfficiencyMap[0] = std::move(formula);
    }
  }

  void Init() override
  {
    fParticleInputArray = ImportArray(Steer<std::string>("ParticleInputArray", "Delphes/allParticles"));
    fPartonInputArray = ImportArray(Steer<std::string>("PartonInputArray", "Delphes/partons"));
    fTrackInputArray = ImportArray(Steer<std::string>("TrackInputArray", "TrackMerger/tracks"));
    fJetInputArray = ImportArray(Steer<std::string>("JetInputArray", "FastJetFinder/jets"));

    fClassifier = std::make_unique<TrackCountingTauTaggingPartonClassifier>(fParticleInputArray);
    fClassifier->fPTMin = Steer<double>("TauPTMin", 1.0);
    fClassifier->fEtaMax = Steer<double>("TauEtaMax", 2.5);

    fFilter = std::make_unique<DelphesFilter>(fPartonInputArray);
  }
  void Process() override;

private:
  const int fBitNumber;

  const double fDeltaR;
  const double fDeltaRTrack;
  const double fTrackPTMin;

  std::map<int, std::unique_ptr<DelphesFormula> > fEfficiencyMap; //!

  CandidatesCollection fParticleInputArray; //!
  CandidatesCollection fPartonInputArray; //!
  CandidatesCollection fTrackInputArray; //!
  CandidatesCollection fJetInputArray; //!

  std::unique_ptr<TrackCountingTauTaggingPartonClassifier> fClassifier; //!
  std::unique_ptr<DelphesFilter> fFilter;
};

//------------------------------------------------------------------------------

int TrackCountingTauTaggingPartonClassifier::GetCategory(TObject *object)
{
  Candidate *tau = static_cast<Candidate *>(object);

  const TLorentzVector &momentum = tau->Momentum;
  int pdgCode, i, j;

  pdgCode = std::abs(tau->PID);
  if(pdgCode != 15) return -1;

  if(momentum.Pt() <= fPTMin || std::fabs(momentum.Eta()) > fEtaMax) return -1;

  if(tau->D1 < 0) return -1;

  if(tau->D2 < tau->D1) return -1;

  if(tau->D1 >= static_cast<int>(fParticleInputArray->size()) || tau->D2 >= static_cast<int>(fParticleInputArray->size()))
  {
    throw runtime_error("tau's daughter index is greater than the ParticleInputArray size");
  }

  for(i = tau->D1; i <= tau->D2; ++i)
  {
    Candidate *daughter1 = static_cast<Candidate *>(fParticleInputArray->at(i));
    pdgCode = std::abs(daughter1->PID);
    if(pdgCode == 11 || pdgCode == 13 || pdgCode == 15)
      return -1;
    else if(pdgCode == 24)
    {
      if(daughter1->D1 < 0) return -1;
      for(j = daughter1->D1; j <= daughter1->D2; ++j)
      {
        Candidate *daughter2 = static_cast<Candidate *>(fParticleInputArray->at(j));
        pdgCode = std::abs(daughter2->PID);
        if(pdgCode == 11 || pdgCode == 13) return -1;
      }
    }
  }

  return 0;
}

//------------------------------------------------------------------------------

void TrackCountingTauTagging::Process()
{
  TLorentzVector tauMomentum;
  double pt, eta, phi, e;
  map<int, std::unique_ptr<DelphesFormula> >::iterator itEfficiencyMap;
  int charge, i, identifier;

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
        if(std::abs(daughter->PID) == 16) continue;
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
