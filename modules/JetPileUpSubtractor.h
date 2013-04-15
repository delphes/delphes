#ifndef JetPileUpSubtractor_h
#define JetPileUpSubtractor_h

/** \class JetPileUpSubtractor
 *
 *  Subtract pile-up contribution from jets using the fastjet area method
 *
 *  $Date: 2012-11-18 15:57:08 +0100 (Sun, 18 Nov 2012) $
 *  $Revision: 814 $
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <deque>

class TObjArray;
class DelphesFormula;

class JetPileUpSubtractor: public DelphesModule
{
public:

  JetPileUpSubtractor();
  ~JetPileUpSubtractor();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fJetPTMin;

  TIterator *fItJetInputArray; //!

  const TObjArray *fJetInputArray; //!
  const TObjArray *fRhoInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(JetPileUpSubtractor, 1)
};

#endif
