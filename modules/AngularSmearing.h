#ifndef AngularSmearing_h
#define AngularSmearing_h

/** \class AngularSmearing
 *
 *  Performs transverse angular resolution smearing.
 *
 *  $Date: 2014-06-17 16:58:53 +0100  $
 *  
 *
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;
class DelphesFormula;

class AngularSmearing: public DelphesModule
{
public:

  AngularSmearing();
  ~AngularSmearing();

  void Init();
  void Process();
  void Finish();

private:

  DelphesFormula *fFormulaEta; //!
  DelphesFormula *fFormulaPhi; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!
  
  TObjArray *fOutputArray; //!

  ClassDef(AngularSmearing, 1)
};

#endif
