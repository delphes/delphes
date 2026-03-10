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
  TIterator *fItInputArray{nullptr}; //!

  const TObjArray *fInputArray{nullptr}; //!

  TObjArray *fOutputArray{nullptr}; //!

  ClassDef(BeamSpotFilter, 1)
};

#endif
