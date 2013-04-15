#ifndef ExRootTreeWriter_h
#define ExRootTreeWriter_h

/** \class ExRootTreeWriter
 *
 *  Class handling output ROOT tree
 *
 *  $Date: 2008-06-04 13:57:27 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "TNamed.h"

#include <set>

class TFile;
class TTree;
class TClass;
class ExRootTreeBranch;

class ExRootTreeWriter : public TNamed
{
public:

  ExRootTreeWriter(TFile *file = 0, const char *treeName = "Analysis");
  ~ExRootTreeWriter();

  void SetTreeFile(TFile *file) { fFile = file; }
  void SetTreeName(const char *name) { fTreeName = name; }

  ExRootTreeBranch *NewBranch(const char *name, TClass *cl);

  void Clear();
  void Fill();
  void Write();

private:

  TTree *NewTree();

  TFile *fFile; //!
  TTree *fTree; //!

  TString fTreeName; //!

  std::set<ExRootTreeBranch*> fBranches; //!

  ClassDef(ExRootTreeWriter, 1)
};

#endif /* ExRootTreeWriter */

