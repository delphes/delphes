#ifndef EnergySmearing_h
#define EnergySmearing_h

/** \class EnergySmearing
 *
 *  Performs energy resolution smearing.
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

class EnergySmearing: public DelphesModule
{
public:

  EnergySmearing();
  ~EnergySmearing();

  void Init();
  void Process();
  void Finish();

private:

  DelphesFormula *fFormula; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!
  
  TObjArray *fOutputArray; //!

  ClassDef(EnergySmearing, 1)
};

#endif
