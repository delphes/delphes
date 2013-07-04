#ifndef PileUpMergerPythia8_h
#define PileUpMergerPythia8_h

/** \class PileUpMergerPythia8
 *
 *  Merges particles from pile-up sample into event
 *
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TObjArray;

namespace Pythia8
{
 class Pythia;
};

class PileUpMergerPythia8: public DelphesModule
{
public:

  PileUpMergerPythia8();
  ~PileUpMergerPythia8();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fMeanPileUp;
  Double_t fZVertexSpread;
  Double_t fPTMin;

  Pythia8::Pythia *fPythia; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(PileUpMergerPythia8, 1)
};

#endif
