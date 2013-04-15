#ifndef Merger_h
#define Merger_h

/** \class Merger
 *
 *  Merges multiple input arrays into one output array
 *  and sums transverse momenta of all input objects.
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <vector>

class TIterator;
class TObjArray;
class DelphesFormula;

class Merger: public DelphesModule
{
public:

  Merger();
  ~Merger();

  void Init();
  void Process();
  void Finish();

private:

  std::vector< TIterator * > fInputList; //!

  TObjArray *fOutputArray; //!
  TObjArray *fMomentumOutputArray; //!
  TObjArray *fEnergyOutputArray; //!

  ClassDef(Merger, 1)
};

#endif
