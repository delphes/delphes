#ifndef ImpactParameterSmearing_h
#define ImpactParameterSmearing_h

/** \class ImpactParameterSmearing
 *
 *  Performs transverse impact parameter smearing.
 *
 *  $Date: 2014-16-03 14:57:44 +0100   
 *
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;
class DelphesFormula;

class ImpactParameterSmearing: public DelphesModule
{
public:

  ImpactParameterSmearing();
  ~ImpactParameterSmearing();

  void Init();
  void Process();
  void Finish();

private:

  DelphesFormula *fFormula; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!
  
  TObjArray *fOutputArray; //!

  ClassDef(ImpactParameterSmearing, 1)
};

#endif
