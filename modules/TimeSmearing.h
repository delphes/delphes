#ifndef TimeSmearing_h
#define TimeSmearing_h

/** \class TimeSmearing
 *
 *  Performs transverse time smearing.
 *
 *  $Date: 2013-12-12 14:57:44 +0100 (Tue, 12 Dec 2013) $
 *
 *  \author Michele Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;
class DelphesFormula;

class TimeSmearing: public DelphesModule
{
public:

  TimeSmearing();
  ~TimeSmearing();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fTimeResolution;
 
  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!
  
  TObjArray *fOutputArray; //!

  ClassDef(TimeSmearing, 1)
};

#endif
