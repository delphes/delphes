
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

#include "classes/DelphesFactory.h"
#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeBranch.h"

#include "TClass.h"
#include "TObjArray.h"

using namespace std;

//------------------------------------------------------------------------------

DelphesFactory::DelphesFactory(const char *name) :
  TNamed(name, ""), fObjArrays(0)
{
  fObjArrays = new ExRootTreeBranch("PermanentObjArrays", TObjArray::Class(), 0);
}

//------------------------------------------------------------------------------

DelphesFactory::~DelphesFactory()
{
  if(fObjArrays) delete fObjArrays;

  map< const TClass*, ExRootTreeBranch* >::iterator itBranches;
  for(itBranches = fBranches.begin(); itBranches != fBranches.end(); ++itBranches)
  {
    delete (itBranches->second);
  }
}

//------------------------------------------------------------------------------

void DelphesFactory::Clear()
{
  set<TObject *>::iterator itPool;
  for(itPool = fPool.begin(); itPool != fPool.end(); ++itPool)
  {
    (*itPool)->Clear();
  }

  TProcessID::SetObjectCount(0);

  map< const TClass*, ExRootTreeBranch* >::iterator itBranches;
  for(itBranches = fBranches.begin(); itBranches != fBranches.end(); ++itBranches)
  {
    itBranches->second->Clear();
  }
}

//------------------------------------------------------------------------------

TObjArray *DelphesFactory::NewPermanentArray()
{
  TObjArray *array = static_cast<TObjArray *>(fObjArrays->NewEntry());
  fPool.insert(array);
  return array;
}

//------------------------------------------------------------------------------

Candidate *DelphesFactory::NewCandidate()
{
  Candidate *object = New<Candidate>();
  object->SetFactory(this);
  TProcessID::AssignID(object);
  return object;
}

//------------------------------------------------------------------------------

TObject *DelphesFactory::New(TClass *cl)
{
  TObject *object = 0;
  ExRootTreeBranch *branch = 0;
  map<const TClass *, ExRootTreeBranch *>::iterator it = fBranches.find(cl);

  if(it != fBranches.end())
  {
    branch = it->second;
  }
  else
  {
    branch = new ExRootTreeBranch(cl->GetName(), cl, 0);
    fBranches.insert(make_pair(cl, branch));
  }

  object = branch->NewEntry();
  object->Clear();
  return object;
}

//------------------------------------------------------------------------------

