/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "classes/DelphesCylindricalFormula.h"

#include "TString.h"

#include <stdexcept>
#include <string>

using namespace std;

//------------------------------------------------------------------------------

DelphesCylindricalFormula::DelphesCylindricalFormula() :
  TFormula()
{
}

//------------------------------------------------------------------------------

DelphesCylindricalFormula::DelphesCylindricalFormula(const char *name, const char *expression) :
  TFormula()
{
}

//------------------------------------------------------------------------------

DelphesCylindricalFormula::~DelphesCylindricalFormula()
{
}

//------------------------------------------------------------------------------

Int_t DelphesCylindricalFormula::Compile(const char *expression)
{
  string buffer;
  const char *it;
  for(it = expression; *it; ++it)
  {
    if(*it == ' ' || *it == '\t' || *it == '\r' || *it == '\n' || *it == '\\' ) continue;
    buffer.push_back(*it);
  }
  if(TFormula::Compile(buffer.c_str()) != 0)
  {
    throw runtime_error("Invalid formula.");
  }
  return 0;
}

//------------------------------------------------------------------------------

Double_t DelphesCylindricalFormula::Eval(Double_t r, Double_t phi, Double_t z)
{
   Double_t x[3] = {r, phi, z};
   return EvalPar(x);
}

//------------------------------------------------------------------------------

Int_t DelphesCylindricalFormula::DefinedVariable(TString &chaine, Int_t &action)
{
  action = kVariable;
  if(chaine == "r")
  {
    if(fNdim < 1) fNdim = 1;
    return 0;
  }
  else if(chaine == "phi")
  {
    if(fNdim < 2) fNdim = 2;
    return 1;
  }
  else if(chaine == "z")
  {
    if(fNdim < 3) fNdim = 3;
    return 2;
  }
  return -1;
}

//------------------------------------------------------------------------------
