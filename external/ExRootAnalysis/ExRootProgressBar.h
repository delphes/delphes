#ifndef ExRootProgressBar_h
#define ExRootProgressBar_h

#include "Rtypes.h"

class ExRootProgressBar
{
public:

  ExRootProgressBar(Long64_t entries, Int_t width = 64);
  ~ExRootProgressBar();

  void Update(Long64_t entry, Long64_t eventCounter = 0, Bool_t last = kFALSE);
  void Finish();

private:

  Long64_t fEntries, fEventCounter;
  Int_t fWidth;

  ULong64_t fTime;
  Int_t fHashes;

  char *fBar; //!

};

#endif /* ExRootProgressBar */

