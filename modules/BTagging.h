#ifndef BTagging_h
#define BTagging_h

/** \class BTagging
 *
 *  Determines origin of jet,
 *  applies b-tagging efficiency (miss identification rate) formulas
 *  and sets b-tagging flags 
 *
 *  $Date$
 *  $Revision$
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
class BTaggingPartonClassifier;

class BTagging: public DelphesModule
{
public:

  BTagging();
  ~BTagging();

  void Init();
  void Process();
  void Finish();

private:

  Int_t fBitNumber;

  Double_t fDeltaR;

  std::map< Int_t, DelphesFormula * > fEfficiencyMap; //!
  
  BTaggingPartonClassifier *fClassifier; //!
  
  ExRootFilter *fFilter;

  TIterator *fItPartonInputArray; //!
  
  TIterator *fItJetInputArray; //!

  const TObjArray *fPartonInputArray; //!
  
  const TObjArray *fJetInputArray; //!

  ClassDef(BTagging, 1)
};

#endif
