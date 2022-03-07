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

#include "classes/DelphesCscClusterFormula.h"
#include "classes/DelphesClasses.h"

#include "TString.h"

#include <stdexcept>
#include <iostream>

using namespace std;

//------------------------------------------------------------------------------

DelphesCscClusterFormula::DelphesCscClusterFormula() :
  TFormula()
{
}

//------------------------------------------------------------------------------

DelphesCscClusterFormula::DelphesCscClusterFormula(const char *name, const char *expression) :
  TFormula()
{
}

//------------------------------------------------------------------------------

DelphesCscClusterFormula::~DelphesCscClusterFormula()
{
}

//------------------------------------------------------------------------------

Int_t DelphesCscClusterFormula::Compile(const char *expression)
{
  TString buffer;
  const char *it;
  for(it = expression; *it; ++it)
  {
    if(*it == ' ' || *it == '\t' || *it == '\r' || *it == '\n' || *it == '\\') continue;
    buffer.Append(*it);
  }
  buffer.ReplaceAll("decayR", "x");
  buffer.ReplaceAll("decayZ", "y");
  buffer.ReplaceAll("Ehad", "z");
  buffer.ReplaceAll("Eem", "t");

#if ROOT_VERSION_CODE < ROOT_VERSION(6, 3, 0)
  TFormula::SetMaxima(100000, 1000, 1000000);
#endif
  if(TFormula::Compile(buffer) != 0)
  {
    throw runtime_error("Invalid formula.");
  }
  return 0;
}

//------------------------------------------------------------------------------

Double_t DelphesCscClusterFormula::Eval(Double_t decayR, Double_t decayZ, Double_t Ehad, Double_t Eem)
{
  Double_t x[4] = {decayR, decayZ, Ehad, Eem};
  return EvalPar(x);
}

//------------------------------------------------------------------------------
