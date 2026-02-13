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

#include "classes/DelphesModel.h"
#include "classes/DelphesModule.h"

class Candidate;

class BeamSpotFilter : public DelphesModule
{
public:
  BeamSpotFilter() = default;

  void Init();
  void Process();
  void Finish();

private:
  Float_t fPassedOne;

  InputHandle<std::vector<Candidate> > fInputArray; //!
  OutputHandle<std::vector<Candidate> > fOutputArray; //!

  ClassDef(BeamSpotFilter, 1)
};

#endif
