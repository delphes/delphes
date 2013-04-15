#ifndef MomentumSmearing_h
#define MomentumSmearing_h

/** \class MomentumSmearing
 *
 *  Performs transverse momentum resolution smearing.
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;
class DelphesFormula;

class MomentumSmearing: public DelphesModule
{
public:

  MomentumSmearing();
  ~MomentumSmearing();

  void Init();
  void Process();
  void Finish();

private:

  DelphesFormula *fFormula; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!
  
  TObjArray *fOutputArray; //!

  ClassDef(MomentumSmearing, 1)
};

#endif
