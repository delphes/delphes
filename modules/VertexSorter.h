#ifndef VertexSorter_h
#define VertexSorter_h

/** \class VertexSorter
 *
 *  Merges particles from pile-up sample into event
 *
 *
 *  $Date: 2013-02-12 15:13:59 +0100 (Tue, 12 Feb 2013) $
 *  $Revision: 907 $
 *
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"
#include "classes/DelphesClasses.h"
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

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

  static bool secondDescending (pair<unsigned, double>, pair<unsigned, double>);
  static bool secondAscending (pair<unsigned, double>, pair<unsigned, double>);

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
