#ifndef ExRootTreeBranch_h
#define ExRootTreeBranch_h

/** \class ExRootTreeBranch
 *
 *  Class handling object creation.
 *  It is also used for output ROOT tree branches
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "Rtypes.h"

class TTree;
class TClonesArray;

class ExRootTreeBranch
{
public:

  ExRootTreeBranch(const char *name, TClass *cl, TTree *tree = 0);
  ~ExRootTreeBranch();

  TObject *NewEntry();
  void Clear();

private:

  Int_t fSize, fCapacity; //!
  TClonesArray *fData; //!
};

#endif /* ExRootTreeBranch */

