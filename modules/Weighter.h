#ifndef Weighter_h
#define Weighter_h

/** \class Weighter
 *
 *  Apply a weight depending on PDG code.
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <set>
#include <map>

class TObjArray;

class Weighter: public DelphesModule
{
public:

  Weighter();
  ~Weighter();

  void Init();
  void Process();
  void Finish();

private:

#if !defined(__CINT__) && !defined(__CLING__)
  struct TIndexStruct
  {
    Int_t codes[4];
    bool operator< (const TIndexStruct &value) const;
  };

  std::set<Int_t> fWeightSet, fCodeSet;
  std::map<TIndexStruct, Double_t> fWeightMap;
#endif

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(Weighter, 1)
};

#endif
