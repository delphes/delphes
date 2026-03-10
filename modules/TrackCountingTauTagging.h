#ifndef TrackCountingTauTagging_h
#define TrackCountingTauTagging_h

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

#include "classes/DelphesModule.h"

#include <map>

class TObjArray;
class DelphesFormula;

class ExRootFilter;
class TrackCountingTauTaggingPartonClassifier;

class TrackCountingTauTagging : public DelphesModule
{
public:
  TrackCountingTauTagging();
  ~TrackCountingTauTagging();

  void Init();
  void Process();
  void Finish();

private:
  Int_t fBitNumber;

  Double_t fDeltaR;
  Double_t fDeltaRTrack;
  Double_t fTrackPTMin;

  std::map<Int_t, std::unique_ptr<DelphesFormula> > fEfficiencyMap; //!

  std::unique_ptr<TrackCountingTauTaggingPartonClassifier> fClassifier; //!
  std::unique_ptr<ExRootFilter> fFilter;

  TIterator *fItPartonInputArray{nullptr}; //!
  TIterator *fItTrackInputArray{nullptr}; //!
  TIterator *fItJetInputArray{nullptr}; //!

  const TObjArray *fParticleInputArray{nullptr}; //!
  const TObjArray *fTrackInputArray{nullptr}; //!
  const TObjArray *fPartonInputArray{nullptr}; //!
  const TObjArray *fJetInputArray{nullptr}; //!

  ClassDef(TrackCountingTauTagging, 1)
};

#endif
