#include "classes/DelphesTF2.h"
#include "TString.h"
#include <stdexcept>
#include <string>

using namespace std;

//------------------------------------------------------------------------------

DelphesTF2::DelphesTF2() :
  TF2()
{
}

//------------------------------------------------------------------------------

DelphesTF2::DelphesTF2(const char *name, const char *expression) :
  TF2(name,expression)
{
}

//------------------------------------------------------------------------------

DelphesTF2::~DelphesTF2()
{
}

//------------------------------------------------------------------------------
/*
Int_t DelphesTF2::Compile(const char *expression)
{
  string buffer;
  const char *it;
  for(it = expression; *it; ++it)
  {
    if(*it == ' ' || *it == '\t' || *it == '\r' || *it == '\n' || *it == '\\' ) continue;
    buffer.push_back(*it);
  }
  
  TFormula f("",buffer.c_str());
  if(f.Compile() != 0)
  {
    throw runtime_error("Invalid formula.");
  }
  return 0;
}

//------------------------------------------------------------------------------

Double_t DelphesTF2::Eval(Double_t z, Double_t t)
{
   Double_t x[2] = {z , t};
   return EvalPar(x);
}

//------------------------------------------------------------------------------

void DelphesTF2::ReplaceAll(std::string& subject, const std::string& search, const std::string& replace)
{
   size_t pos = 0;
   while ((pos = subject.find(search, pos)) != std::string::npos) 
   {
      subject.replace(pos, search.length(), replace);
      pos += replace.length();
   }
}


//------------------------------------------------------------------------------

string DelphesTF2::ChangeVariables(const char *expression)
{
   DelphesTF2::Compile(expression);
   
   string str = string(expression);
   ReplaceAll(str,"z","x");
   ReplaceAll(str,"t","y");
   
   return str;
}

*/
//------------------------------------------------------------------------------

Int_t DelphesTF2::DefinedVariable(TString &chaine, Int_t &action)
{
  action = kVariable;
  if(chaine == "z")
  {
    if(fNdim < 1) fNdim = 1;
    return 0;
  }
  else if(chaine == "t")
  {
    if(fNdim < 2) fNdim = 2;
    return 1;
  }
  return -1;
}

//------------------------------------------------------------------------------
