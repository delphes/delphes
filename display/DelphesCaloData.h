#ifndef DelphesCaloData_h
#define DelphesCaloData_h

#include "TEveCaloData.h"

class DelphesCaloData: public TEveCaloDataVec 
{
public:
  
  DelphesCaloData(Int_t nslices);

  ~DelphesCaloData();

  void ClearTowers();

  ClassDef(DelphesCaloData, 1)
};

#endif /* DelphesCaloData_h */

