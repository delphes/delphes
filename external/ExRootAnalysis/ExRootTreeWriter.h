#ifndef ExRootTreeWriter_h
#define ExRootTreeWriter_h

/** \class ExRootTreeWriter
 *
 *  Class handling output ROOT tree
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

  virtual void Clear(Option_t *option = "");
  void Fill();
  virtual Int_t Write(const char *name=0, Int_t option=0, Int_t bufsize=0) const;
  virtual Int_t Write(const char *name=0, Int_t option=0, Int_t bufsize=0);

private:

  TTree *NewTree();

  TFile *fFile; //!
  TTree *fTree; //!

  TString fTreeName; //!

  std::set<ExRootTreeBranch*> fBranches; //!

  ClassDef(ExRootTreeWriter, 1)
};

#endif /* ExRootTreeWriter */

