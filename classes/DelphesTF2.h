#ifndef DelphesTF2_h
#define DelphesTF2_h

#include "TF2.h"
#include "TFormula.h"

#include <string>

class DelphesTF2: public TF2
{
public:

  DelphesTF2();

  DelphesTF2(const char *name, const char *expression);

  ~DelphesTF2();

//  Int_t Compile(const char *expression);

//  Double_t Eval(Double_t z, Double_t t);
  
//  std::string ChangeVariables(const char *expression);
  
//  void ReplaceAll(std::string& str, const std::string& from, const std::string& to);
  
  Int_t DefinedVariable(TString &variable, Int_t &action);

};

#endif /* DelphesTF2_h */

