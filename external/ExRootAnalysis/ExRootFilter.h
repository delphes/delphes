#ifndef ExRootFilter_h
#define ExRootFilter_h

#include "Rtypes.h"

#include <map>

class ExRootClassifier;
class TSeqCollection;
class TObjArray;
class TIterator;

class ExRootFilter
{
public:

  ExRootFilter(const TSeqCollection *collection);
  ~ExRootFilter();

  void Reset(ExRootClassifier *classifier = 0);

  TObjArray *GetSubArray(ExRootClassifier *classifier, Int_t category);

private:

  const TSeqCollection *fCollection; //!
  TIterator *fIter; //!

  std::map<ExRootClassifier*, std::pair<Bool_t, std::map<Int_t, TObjArray*> > > fMap; //!

};

#endif /* ExRootFilter */

