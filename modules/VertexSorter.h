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
  TObjArray *fInputArray{nullptr};

  TObjArray *fTrackInputArray{nullptr};
  TObjArray *fJetInputArray{nullptr};
  TObjArray *fBeamSpotInputArray{nullptr};

  std::unique_ptr<TIterator> fItTrackInputArray;
  std::unique_ptr<TIterator> fItJetInputArray;

  TObjArray *fOutputArray{nullptr};

  std::string fMethod;

  ClassDef(VertexSorter, 1)
};

#endif
