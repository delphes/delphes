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
  TObjArray *fInputArray = nullptr;

  TObjArray *fTrackInputArray = nullptr;
  TIterator *fItTrackInputArray = nullptr;

  TObjArray *fJetInputArray = nullptr;
  TIterator *fItJetInputArray = nullptr;

  TObjArray *fBeamSpotInputArray = nullptr;
  TIterator *fItBeamSpotInputArray = nullptr;

  TObjArray *fOutputArray = nullptr;

  std::string fMethod;

  ClassDef(VertexSorter, 1)
};

#endif
