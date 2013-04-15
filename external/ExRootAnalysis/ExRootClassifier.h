#ifndef ExRootClassifier_h
#define ExRootClassifier_h

#include <Rtypes.h>

class TObject;

class ExRootClassifier
{
public:
  virtual ~ExRootClassifier() {}
  virtual Int_t GetCategory(TObject *object) = 0;

};

#endif /* ExRootClassifier */

