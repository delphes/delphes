#include "classes/DelphesTF2.h"
#include "TString.h"
#include <stdexcept>
#include <string>

using namespace std;

//------------------------------------------------------------------------------

DelphesTF2::DelphesTF2() :
  TF2()
{
}

//------------------------------------------------------------------------------

DelphesTF2::DelphesTF2(const char *name, const char *expression) :
  TF2(name,expression)
{
}

//------------------------------------------------------------------------------

DelphesTF2::~DelphesTF2()
{
}

//------------------------------------------------------------------------------

Int_t DelphesTF2::DefinedVariable(TString &chaine, Int_t &action)
{
  action = kVariable;
  if(chaine == "z")
  {
    if(fNdim < 1) fNdim = 1;
    return 0;
  }
  else if(chaine == "t")
  {
    if(fNdim < 2) fNdim = 2;
    return 1;
  }
  return -1;
}

//------------------------------------------------------------------------------
