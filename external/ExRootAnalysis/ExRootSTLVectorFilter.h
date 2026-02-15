#ifndef ExRootFilter_h
#define ExRootFilter_h

#include "Rtypes.h"

#include <map>

class Candidate;
class ExRootClassifier;

class ExRootSTLVectorFilter
{
public:
  ExRootSTLVectorFilter(const std::vector<Candidate> &collection);

  void Reset(ExRootClassifier *classifier = 0);

  std::vector<Candidate> GetSubArray(ExRootClassifier *classifier, Int_t category);

private:
  const std::vector<Candidate> &fCollection; //!

  using TCategoryMap = std::map<Int_t, std::vector<Candidate> >;
  using TClassifierMap = std::map<ExRootClassifier *, std::pair<Bool_t, TCategoryMap> >;

  TClassifierMap fMap; //!
};

#endif /* ExRootFilter */
