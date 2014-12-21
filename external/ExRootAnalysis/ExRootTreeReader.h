#ifndef ExRootTreeReader_h
#define ExRootTreeReader_h

/** \class ExRootTreeReader
 *
 *  Class simplifying access to ROOT tree branches
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "TROOT.h"
#include "TNamed.h"
#include "TChain.h"
#include "TFile.h"

#include <map>

class ExRootTreeReader : public TNamed
{
public :

  ExRootTreeReader(TTree *tree = 0);
  ~ExRootTreeReader();

  void SetTree(TTree *tree) { fChain = tree; }

  Long64_t GetEntries() const { return fChain ? static_cast<Long64_t>(fChain->GetEntries()) : 0; }
  Bool_t ReadEntry(Long64_t entry);

  TClonesArray *UseBranch(const char *branchName);

private:

  Bool_t Notify();

  TTree *fChain; //! pointer to the analyzed TTree or TChain
  Int_t fCurrentTree; //! current Tree number in a TChain

  typedef std::map<TString, std::pair<TBranch*, TClonesArray*> > TBranchMap;

  TBranchMap fBranchMap; //!

  ClassDef(ExRootTreeReader, 1)
};

#endif // ExRootTreeReader_h
