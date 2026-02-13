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

#include "classes/DelphesModel.h"
#include "classes/DelphesModule.h"

#include <map>

class Candidate;
class DelphesFormula;

class ExRootSTLVectorFilter;
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

  std::map<Int_t, DelphesFormula *> fEfficiencyMap; //!

  TrackCountingTauTaggingPartonClassifier *fClassifier; //!

  ExRootSTLVectorFilter *fFilter;

  InputHandle<std::vector<Candidate> > fParticleInputArray; //!
  InputHandle<std::vector<Candidate> > fTrackInputArray; //!
  InputHandle<std::vector<Candidate> > fPartonInputArray; //!
  InputHandle<std::vector<Candidate> > fJetInputArray; //!

  ClassDef(TrackCountingTauTagging, 1)
};

#endif
