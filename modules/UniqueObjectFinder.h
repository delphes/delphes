#ifndef UniqueObjectFinder_h
#define UniqueObjectFinder_h

/** \class UniqueObjectFinder
 *
 *  Finds uniquely identified photons, electrons, taus and jets.
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <map>

class TIterator;
class TObjArray;
class Candidate;

class UniqueObjectFinder: public DelphesModule
{
public:

  UniqueObjectFinder();
  ~UniqueObjectFinder();

  void Init();
  void Process();
  void Finish();

private:

  Bool_t Unique(Candidate *candidate, std::map< TIterator *, TObjArray * >::iterator itInputMap);

  std::map< TIterator *, TObjArray * > fInputMap; //!

  ClassDef(UniqueObjectFinder, 1)
};

#endif
