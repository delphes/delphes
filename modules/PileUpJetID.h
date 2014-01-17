#ifndef PileUpJetID_h
#define PileUpJetID_h

/** \class PileUpJetID
 *
 *  CMS PileUp Jet ID Variables, based on http://cds.cern.ch/record/1581583
 *
 *  \author S. Zenz, December 2013
 *  
 *
 */


#include "classes/DelphesModule.h"

#include <deque>

class TObjArray;
class DelphesFormula;

class PileUpJetID: public DelphesModule
{
public:

  PileUpJetID();
  ~PileUpJetID();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fJetPTMin;
  Double_t fParameterR;

  // If set to true, may have weird results for PFCHS
  // If set to false, uses everything within dR < fParameterR even if in other jets &c.
  // Results should be very similar for PF
  Int_t fUseConstituents; 

  Bool_t fAverageEachTower;

  TIterator *fItJetInputArray; //!

  const TObjArray *fJetInputArray; //!

  const TObjArray *fTrackInputArray; // SCZ
  const TObjArray *fNeutralInputArray; 

  TIterator *fItTrackInputArray; // SCZ
  TIterator *fItNeutralInputArray; // SCZ

  TObjArray *fOutputArray; //!
  
  TIterator *fItVertexInputArray; //!
  const TObjArray *fVertexInputArray; //!

  Double_t fZVertexResolution;

  ClassDef(PileUpJetID, 1)
};

#endif

