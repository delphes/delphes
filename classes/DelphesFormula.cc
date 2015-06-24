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


#include "classes/DelphesFormula.h"

#include "TString.h"

#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

DelphesFormula::DelphesFormula() :
  TFormula()
{
}

//------------------------------------------------------------------------------

DelphesFormula::DelphesFormula(const char *name, const char *expression) :
  TFormula()
{
}

//------------------------------------------------------------------------------

DelphesFormula::~DelphesFormula()
{
}

//------------------------------------------------------------------------------

Int_t DelphesFormula::Compile(const char *expression)
{
  TString buffer;
  const char *it;
  for(it = expression; *it; ++it)
  {
    if(*it == ' ' || *it == '\t' || *it == '\r' || *it == '\n' || *it == '\\' ) continue;
    buffer.Append(*it);
  }
  buffer.ReplaceAll("pt", "x");
  buffer.ReplaceAll("eta", "y");
  buffer.ReplaceAll("phi", "z");
  buffer.ReplaceAll("energy", "t");
  if(TFormula::Compile(buffer) != 0)
  {
    throw runtime_error("Invalid formula.");
  }
  return 0;
}

//------------------------------------------------------------------------------

Double_t DelphesFormula::Eval(Double_t pt, Double_t eta, Double_t phi, Double_t energy)
{
   Double_t x[4] = {pt, eta, phi, energy};
   return EvalPar(x);
}

//------------------------------------------------------------------------------
