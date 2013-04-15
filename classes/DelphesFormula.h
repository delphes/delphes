#ifndef DelphesFormula_h
#define DelphesFormula_h

#include "TFormula.h"

class DelphesFormula: public TFormula
{
public:

  DelphesFormula();

  DelphesFormula(const char *name, const char *expression);

  ~DelphesFormula();

  Int_t Compile(const char *expression);

  Double_t Eval(Double_t pt, Double_t eta = 0, Double_t phi = 0, Double_t energy = 0);

  Int_t DefinedVariable(TString &variable, Int_t &action);
};

#endif /* DelphesFormula_h */

