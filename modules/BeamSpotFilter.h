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
  const TObjArray *fInputArray{nullptr}; //!
  std::unique_ptr<TIterator> fItInputArray; //!

  TObjArray *fOutputArray{nullptr}; //!

  ClassDef(BeamSpotFilter, 1)
};

#endif
