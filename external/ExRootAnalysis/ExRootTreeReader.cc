
/** \class ExRootTreeReader
 *
 *  Class simplifying access to ROOT tree branches
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootTreeReader.h"

#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TBranchElement.h"

#include <iostream>

using namespace std;

//------------------------------------------------------------------------------

ExRootTreeReader::ExRootTreeReader(TTree *tree) :
  fChain(tree), fCurrentTree(-1)
{
}

//------------------------------------------------------------------------------

ExRootTreeReader::~ExRootTreeReader()
{
  TBranchMap::iterator itBranchMap;

  for(itBranchMap = fBranchMap.begin(); itBranchMap != fBranchMap.end(); ++itBranchMap)
  {
    delete itBranchMap->second.second;
  }
}

//------------------------------------------------------------------------------

Bool_t ExRootTreeReader::ReadEntry(Long64_t entry)
{
  // Read contents of entry.
  if(!fChain) return kFALSE;

  Int_t treeEntry = fChain->LoadTree(entry);
  if(treeEntry < 0) return kFALSE;

  if(fChain->IsA() == TChain::Class())
  {
    TChain *chain = static_cast<TChain*>(fChain);
    if(chain->GetTreeNumber() != fCurrentTree)
    {
      fCurrentTree = chain->GetTreeNumber();
      Notify();
    }
  }

  TBranchMap::iterator itBranchMap;
  TBranch *branch;

  for(itBranchMap = fBranchMap.begin(); itBranchMap != fBranchMap.end(); ++itBranchMap)
  {
    branch = itBranchMap->second.first;
    if(branch)
    {
      branch->GetEntry(treeEntry);
    }
  }

  return kTRUE;
}

//------------------------------------------------------------------------------

TClonesArray *ExRootTreeReader::UseBranch(const char *branchName)
{
  TClonesArray *array = 0;

  TBranchMap::iterator itBranchMap = fBranchMap.find(branchName);

  if(itBranchMap != fBranchMap.end())
  {
    cout << "** WARNING: branch '" << branchName << "' is already in use" << endl;
    array = itBranchMap->second.second;
  }
  else
  {
    TBranch *branch = fChain->GetBranch(branchName);
    if(branch)
    {
      if(branch->IsA() == TBranchElement::Class())
      {
        TBranchElement *element = static_cast<TBranchElement*>(branch);
        const char *className = element->GetClonesName();
        Int_t size = element->GetMaximum();
        TClass *cl = gROOT->GetClass(className);
        if(cl)
        {
          array = new TClonesArray(cl, size);
          array->SetName(branchName);
          fBranchMap.insert(make_pair(branchName, make_pair(branch, array)));
          branch->SetAddress(&array);
        }
      }
    }
  }

  if(!array)
  {
    cout << "** WARNING: cannot access branch '" << branchName << "', return NULL pointer" << endl;
  }

  return array;
}

//------------------------------------------------------------------------------

Bool_t ExRootTreeReader::Notify()
{
  // Called when loading a new file.
  // Get branch pointers.
  if(!fChain) return kFALSE;

  TBranchMap::iterator itBranchMap;
  TBranch *branch;

  for(itBranchMap = fBranchMap.begin(); itBranchMap != fBranchMap.end(); ++itBranchMap)
  {
    branch = fChain->GetBranch(itBranchMap->first);
    if(branch)
    {
      itBranchMap->second.first = branch;
      branch->SetAddress(&(itBranchMap->second.second));
    }
    else
    {
      cout << "** WARNING: cannot get branch '" << itBranchMap->first << "'" << endl;
    }
  }
  return kTRUE;
}

//------------------------------------------------------------------------------

