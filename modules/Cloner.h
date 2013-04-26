#ifndef Clone_h
#define Clone_h

/** \class Clone
 *
 *  Clone candidate array
 *
 *  $Date$
 *  $Revision$
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <deque>

class TObjArray;
class DelphesFormula;

class Cloner: public DelphesModule
{
public:

  Cloner();
  ~Cloner();

  void Init();
  void Process();
  void Finish();

private:

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!
  TObjArray *fOutputArray; //!

  ClassDef(Cloner, 1)
};

#endif
