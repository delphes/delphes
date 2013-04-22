
#include "display/DelphesCaloData.h"

//------------------------------------------------------------------------------

DelphesCaloData::DelphesCaloData(Int_t nslices) :
  TEveCaloDataVec(nslices)
{
}

//------------------------------------------------------------------------------

DelphesCaloData::~DelphesCaloData()
{
}

//------------------------------------------------------------------------------

void DelphesCaloData::ClearTowers()
{
   fGeomVec.clear();
}

//------------------------------------------------------------------------------
