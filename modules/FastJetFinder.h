#ifndef FastJetFinder_h
#define FastJetFinder_h

/** \class FastJetFinder
 *
 *  Finds jets using FastJet library.
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
class TIterator;

namespace fastjet {
  class JetDefinition;
  class AreaDefinition;
  class Selector;
}

class FastJetFinder: public DelphesModule
{
public:

  FastJetFinder();
  ~FastJetFinder();

  void Init();
  void Process();
  void Finish();

private:

  void *fPlugin; //!
  void *fRecomb; //!
  void *fNjettinessPlugin; //!
  
  fastjet::JetDefinition *fDefinition; //!

  Int_t fJetAlgorithm;
  Double_t fParameterR;
  Double_t fJetPTMin;
  Double_t fConeRadius;
  Double_t fSeedThreshold;
  Double_t fConeAreaFraction;
  Int_t fMaxIterations;
  Int_t fMaxPairSize;
  Int_t fIratch;
  Int_t fAdjacencyCut;
  Double_t fOverlapThreshold;

  //-- N (sub)jettiness parameters --
  
  Bool_t fComputeNsubjettiness;
  Double_t fBeta;
  Int_t fAxisMode;
  Double_t fRcutOff; 
  Int_t fN ;
  
  // --- FastJet Area method --------

  fastjet::AreaDefinition *fAreaDefinition;
  Int_t fAreaAlgorithm;
  Bool_t  fComputeRho;

  // -- ghost based areas --
  Double_t fGhostEtaMax;
  Int_t fRepeat;
  Double_t fGhostArea;
  Double_t fGridScatter;
  Double_t fPtScatter;
  Double_t fMeanGhostPt;

  // -- voronoi areas --
  Double_t fEffectiveRfact;

  std::map< Double_t, Double_t > fEtaRangeMap; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!
  TObjArray *fRhoOutputArray; //!

  ClassDef(FastJetFinder, 1)
};

#endif
