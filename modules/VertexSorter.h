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

class Candidate;

class VertexSorter : public DelphesModule
{
public:
  VertexSorter() = default;

  void Init();
  void Process();
  void Finish();

private:
  InputHandle<std::vector<Candidate> > fInputArray; //!
  InputHandle<std::vector<Candidate> > fTrackInputArray; //!
  InputHandle<std::vector<Candidate> > fJetInputArray; //!
  InputHandle<std::vector<Candidate> > fBeamSpotInputArray; //!
  OutputHandle<std::vector<Candidate> > fOutputArray; //!

  std::string fMethod;

  ClassDef(VertexSorter, 1)
};

#endif
