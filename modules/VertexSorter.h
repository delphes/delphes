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
#include "classes/DelphesClasses.h"
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>


using namespace std;

class TObjArray;
class Candidate;
class TVector3;

class VertexSorter: public DelphesModule
{
public:

  VertexSorter();
  ~VertexSorter();

  void Init();
  void Process();
  void Finish();

  static Bool_t secondDescending (pair<UInt_t, Double_t>, pair<UInt_t, Double_t>);
  static Bool_t secondAscending (pair<UInt_t, Double_t>, pair<UInt_t, Double_t>);

private:

  TObjArray *fInputArray;

  TObjArray *fTrackInputArray;
  TIterator *fItTrackInputArray;

  TObjArray *fJetInputArray;
  TIterator *fItJetInputArray;

  TObjArray *fBeamSpotInputArray;
  TIterator *fItBeamSpotInputArray;

  TObjArray *fOutputArray;

  string fMethod;

  ClassDef(VertexSorter, 1)
};

#endif
