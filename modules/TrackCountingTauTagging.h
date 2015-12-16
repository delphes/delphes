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

class TrackCountingTauTagging: public DelphesModule
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

  std::map< Int_t, DelphesFormula * > fEfficiencyMap; //!
  
  TrackCountingTauTaggingPartonClassifier *fClassifier; //!
  
  ExRootFilter *fFilter;

  TIterator *fItPartonInputArray; //!
  
  TIterator *fItTrackInputArray; //!
  
  TIterator *fItJetInputArray; //!

  const TObjArray *fParticleInputArray; //!

  const TObjArray *fTrackInputArray; //!

  const TObjArray *fPartonInputArray; //!
  
  const TObjArray *fJetInputArray; //!

  ClassDef(TrackCountingTauTagging, 1)
};

#endif
