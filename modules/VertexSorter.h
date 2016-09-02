#ifndef VertexSorter_h
#define VertexSorter_h

/** \class VertexSorter
 *
 *
 *  Sorts vertices according to different criteria
 *
 *  \authors A. Hart, M. Selvaggi
 *
 *
*/

#include "classes/DelphesModule.h"

#include <string>

class TObjArray;
class TIterator;
class Candidate;

class VertexSorter: public DelphesModule
{
public:

  VertexSorter();
  ~VertexSorter();

  void Init();
  void Process();
  void Finish();

private:

  TObjArray *fInputArray;

  TObjArray *fTrackInputArray;
  TIterator *fItTrackInputArray;

  TObjArray *fJetInputArray;
  TIterator *fItJetInputArray;

  TObjArray *fBeamSpotInputArray;
  TIterator *fItBeamSpotInputArray;

  TObjArray *fOutputArray;

  std::string fMethod;

  ClassDef(VertexSorter, 1)
};

#endif
