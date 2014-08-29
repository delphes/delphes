#ifndef IdentificationMap_h
#define IdentificationMap_h


/** \class IdentificationMap
 *
 *  Converts particles with some PDG code into another particle, according to parametrized probability as fction of pt eta
 given by user. 
 *
 *  $Date: 2014-08-07 14:57:44 +0100 (Thu, 07 Aug 2014) $
 *  $Revision: 905 $
 *
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;
class DelphesFormula;

class IdentificationMap: public DelphesModule
{
public:

  IdentificationMap();
  ~IdentificationMap();

  void Init();
  void Process();
  void Finish();

private:
  
  typedef std::multimap< Int_t, std::pair<Int_t , DelphesFormula * > > TMisIDMap; //!
  
  TMisIDMap fEfficiencyMap;
  
  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(IdentificationMap, 1)
};

#endif
