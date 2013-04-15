#ifndef ConstituentFilter_h
#define ConstituentFilter_h

/** \class ConstituentFilter
 *
 *  Drops all input objects that are not constituents of any jet.
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
#include <map>

class TIterator;
class TObjArray;

class ConstituentFilter: public DelphesModule
{
public:

  ConstituentFilter();
  ~ConstituentFilter();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fJetPTMin;

  std::vector< TIterator * > fInputList; //!

  std::map< TIterator *, TObjArray * > fInputMap; //!

  TObjArray *fOutputArray; //!

  ClassDef(ConstituentFilter, 1)
};

#endif
