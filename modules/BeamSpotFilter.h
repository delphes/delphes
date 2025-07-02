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

  void ProduceVertices();

  Float_t fPassedOne;

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  Int_t fMode;
  Double_t fResolution; // used internally for fMode==1

  TObjArray *fVertexArray; // used internally for fMode==1
  TIterator *fItVertexArray; // used internally for fMode==1

  ClassDef(BeamSpotFilter, 1)
};

#endif
