#ifndef PileUpJetID_h
#define PileUpJetID_h

/** \class PileUpJetID
 *
 *  CMS PileUp Jet ID Variables
 *
 *  \author S. Zenz
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

  Double_t fMeanSqDeltaRMaxBarrel; // |eta| < 1.5
  Double_t fBetaMinBarrel; // |eta| < 2.5
  Double_t fMeanSqDeltaRMaxEndcap; // 1.5 < |eta| < 4.0
  Double_t fBetaMinEndcap; // 1.5 < |eta| < 4.0
  Double_t fMeanSqDeltaRMaxForward; // |eta| > 4.0

  Double_t fNeutralPTMin;
  Double_t fJetPTMinForNeutrals;

  /*
JAY
---

|Eta|<1.5

meanSqDeltaR betaStar SigEff BgdEff
0.13 0.92 96% 8%
0.13 0.95 97% 16%
0.13 0.97 98% 27%

|Eta|>1.5

meanSqDeltaR betaStar SigEff BgdEff
0.14 0.91 95% 15%
0.14 0.94 97% 19%
0.14 0.97 98% 29%

BRYAN
-----

Barrel (MeanSqDR, Beta, sig eff, bg eff): 
0.10, 0.08, 90%, 8%
0.11, 0.12, 90%, 6%
0.13, 0.16, 89%, 5%

Endcap (MeanSqDR, Beta, sig eff, bg eff):
0.07, 0.06, 89%, 4%
0.08, 0.08, 92%, 6%
0.09, 0.08, 95%, 10%
0.10, 0.08, 97%, 13%

SETH GUESSES FOR |eta| > 4.0
----------------------------

MeanSqDeltaR 
0.07 
0.10 
0.14 
0.2
  */

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
  TObjArray *fNeutralsInPassingJets; // SCZ


  ClassDef(PileUpJetID, 2)
};

#endif
