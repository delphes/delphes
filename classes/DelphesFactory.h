#ifndef DelphesFactory_h
#define DelphesFactory_h

/** \class DelphesFactory
 *
 *  Class handling creation of Candidate,
 *  TObjArray and all other objects.
 *
 *  $Date: 2008-06-04 13:57:25 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "TNamed.h"

#include <map>
#include <set>

class TObjArray;
class Candidate;

class ExRootTreeBranch;

class DelphesFactory: public TNamed
{
public:
  
  DelphesFactory(const char *name = "ObjectFactory");
  ~DelphesFactory();

  void Clear();
 
  TObjArray *NewPermanentArray();

  TObjArray *NewArray() { return New<TObjArray>(); }

  Candidate *NewCandidate();

  TObject *New(TClass *cl);

  template<typename T>
  T *New() { return static_cast<T *>(New(T::Class())); }

private:

  ExRootTreeBranch *fObjArrays; //!

  std::map< const TClass*, ExRootTreeBranch* > fBranches; //!
  std::set< TObject* > fPool; //!
  
  ClassDef(DelphesFactory, 1)
};

#endif /* DelphesFactory */

