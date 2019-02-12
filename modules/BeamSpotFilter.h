//------------------------------------------------------------------------------

#ifndef BeamSpotFilter_h
#define BeamSpotFilter_h

/** \class BeamSpotFilter
 *
 *  Extracts beam spot 
 *
 *  \author Michele Selvaggi
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;

class BeamSpotFilter: public DelphesModule
{
public:
  BeamSpotFilter();
  ~BeamSpotFilter();

  void Init();
  void Process();
  void Finish();

private:
  Float_t fPassedOne;

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(BeamSpotFilter, 1)
};

#endif
