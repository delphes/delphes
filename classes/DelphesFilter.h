#ifndef DelphesFilter_h
#define DelphesFilter_h

#include "Rtypes.h"

#include <map>

#include "classes/DelphesClasses.h"

class ExRootClassifier;

class DelphesFilter
{
public:
  DelphesFilter(const CandidatesCollection &collection);
  ~DelphesFilter() = default;

  void Reset(ExRootClassifier *classifier = 0);

  std::vector<Candidate *> GetSubArray(ExRootClassifier *classifier, Int_t category);

private:
  const CandidatesCollection &fCollection; //!

  using TCategoryMap = std::map<Int_t, std::vector<Candidate *> >;
  using TClassifierMap = std::map<ExRootClassifier *, std::pair<Bool_t, TCategoryMap> >;

  TClassifierMap fMap; //!
};

#endif /* DelphesFilter */
