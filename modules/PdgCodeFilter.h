//------------------------------------------------------------------------------

#ifndef PidFilter_h
#define PidFilter_h

/** \class Efficiency
 *
 *  Removes particles with specific pdg codes 
  *
 *  \author M. Selvaggi
 *
 */

#include "classes/DelphesModule.h"
#include <vector>

class TIterator;
class TObjArray;

class PdgCodeFilter: public DelphesModule
{
public:

  PdgCodeFilter();
  ~PdgCodeFilter();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fPTMin; //!
  
  std::vector<Int_t> fPdgCodes; 

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(PdgCodeFilter, 1)
};

#endif
