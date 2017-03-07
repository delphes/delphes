
/** \class ExRootTreeWriter
 *
 *  Class handling output ROOT tree
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include <iostream>
#include <stdexcept>
#include <sstream>

using namespace std;

ExRootTreeWriter::ExRootTreeWriter(TFile *file, const char *treeName) :
  fFile(file), fTree(0), fTreeName(treeName)
{
}

//------------------------------------------------------------------------------

ExRootTreeWriter::~ExRootTreeWriter()
{
  set<ExRootTreeBranch*>::iterator itBranches;
  for(itBranches = fBranches.begin(); itBranches != fBranches.end(); ++itBranches)
  {
    delete (*itBranches);
  }

  if(fTree) delete fTree;
}

//------------------------------------------------------------------------------

ExRootTreeBranch *ExRootTreeWriter::NewBranch(const char *name, TClass *cl)
{
  if(!fTree) fTree = NewTree();
  ExRootTreeBranch *branch = new ExRootTreeBranch(name, cl, fTree);
  fBranches.insert(branch);
  return branch;
}

//------------------------------------------------------------------------------

void ExRootTreeWriter::Fill()
{
  if(fTree) fTree->Fill();
}

//------------------------------------------------------------------------------

Int_t ExRootTreeWriter::Write(const char *name, Int_t option, Int_t bufsize) const
{
  if (fTree != 0) {
    return fTree->GetCurrentFile()->Write(name, option, bufsize);
  }
  return 0;
}

//------------------------------------------------------------------------------

Int_t ExRootTreeWriter::Write(const char *name, Int_t option, Int_t bufsize)
{
  // forward call to const version of this function
  return static_cast<const ExRootTreeWriter*>(this)->Write(name, option, bufsize);
}

//------------------------------------------------------------------------------

void ExRootTreeWriter::Clear(Option_t *option)
{
  set<ExRootTreeBranch*>::iterator itBranches;
  for(itBranches = fBranches.begin(); itBranches != fBranches.end(); ++itBranches)
  {
    (*itBranches)->Clear(option);
  }
}

//------------------------------------------------------------------------------

TTree *ExRootTreeWriter::NewTree()
{
  if(!fFile) return 0;

  TTree *tree = 0;
  TDirectory *dir = gDirectory;

  fFile->cd();
  tree = new TTree(fTreeName, "Analysis tree");
  dir->cd();

  if(!tree)
  {
    throw runtime_error("can't create output ROOT tree");
  }

  tree->SetDirectory(fFile);
  tree->SetAutoSave(10000000);  // autosave when 10 MB written

  return tree;
}
