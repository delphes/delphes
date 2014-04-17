#ifndef TrackCountingBTagging_h
#define TrackCountingBTagging_h

/** \class TrackCountingBTagging
 *
 *  b-tagging algorithm based on counting tracks with large impact parameter
 *
 *  $Date: 2014-03-27 12:39:14 +0200 (Fri, 27 March 2014) $
 *  $Revision: 1099 $
 *
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <map>

class TObjArray;
class DelphesFormula;

class TrackCountingBTagging: public DelphesModule
{
public:

  TrackCountingBTagging();
  ~TrackCountingBTagging();

  void Init();
  void Process();
  void Finish();

private:

  Int_t fBitNumber;

  Double_t fPtMin;
  Double_t fDeltaR;
  Double_t fIPmax;
  Double_t fSigMin;
  Int_t    fNtracks;

  TIterator *fItTrackInputArray; //!
  TIterator *fItJetInputArray; //!

  const TObjArray *fTrackInputArray; //!
  const TObjArray *fJetInputArray; //!

  ClassDef(TrackCountingBTagging, 1)
};

#endif
