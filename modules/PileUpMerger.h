#ifndef PileUpMerger_h
#define PileUpMerger_h

/** \class PileUpMerger
 *
 *  Merges particles from pile-up sample into event
 *
 *
 *  $Date: 2013-02-12 15:13:59 +0100 (Tue, 12 Feb 2013) $
 *  $Revision: 907 $
 *
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TObjArray;
class DelphesPileUpReader;

class PileUpMerger: public DelphesModule
{
public:

  PileUpMerger();
  ~PileUpMerger();

  void Init();
  void Process();
  void Finish();

private:

  Int_t fPileUpDistribution;
  Double_t fMeanPileUp;
  Double_t fZVertexSpread;

  DelphesPileUpReader *fReader;

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(PileUpMerger, 1)
};

#endif
