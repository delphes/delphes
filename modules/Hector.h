#ifndef Hector_h
#define Hector_h

/** \class Hector
 *
 *  Propagates candidates using Hector library.
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;
class H_BeamLine;

class Hector: public DelphesModule
{
public:

  Hector();
  ~Hector();

  void Init();
  void Process();
  void Finish();

private:

  Int_t fDirection;
  
  Double_t fBeamLineLength, fDistance;
  Double_t fSigmaE, fSigmaX, fSigmaY;
  Double_t fEtaMin;
  
  H_BeamLine *fBeamLine;

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(Hector, 1)
};

#endif
