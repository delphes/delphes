#ifndef LeptonDressing_h
#define LeptonDressing_h

/** \class LeptonDressing
 *
 *
 *  \author P. Demin && A. Mertens - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;

class LeptonDressing: public DelphesModule
{
public:

  LeptonDressing();
  ~LeptonDressing();
  
  void Init();
  void Process();
  void Finish();

private:

  Double_t fDeltaR;
  
  TIterator *fItDressingInputArray; //!
  
  TIterator *fItCandidateInputArray; //!

  const TObjArray *fDressingInputArray; //!
  
  const TObjArray *fCandidateInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(LeptonDressing, 1)
};

#endif
