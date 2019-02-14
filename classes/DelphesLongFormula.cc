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

#include "classes/DelphesLongFormula.h"

#include "TString.h"

#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

DelphesLongFormula::DelphesLongFormula() :
  TFormula()
{
}

//------------------------------------------------------------------------------

DelphesLongFormula::DelphesLongFormula(const char *name, const char *expression) :
  TFormula()
{
}

//------------------------------------------------------------------------------

DelphesLongFormula::~DelphesLongFormula()
{
}

//------------------------------------------------------------------------------

Int_t DelphesLongFormula::Compile(const char *expression)
{
  TString buffer;
  const char *it;
  for(it = expression; *it; ++it)
  {
    if(*it == ' ' || *it == '\t' || *it == '\r' || *it == '\n' || *it == '\\') continue;
    buffer.Append(*it);
  }

#if ROOT_VERSION_CODE < ROOT_VERSION(6, 0, 0)
  TFormula::SetMaxima(100000, 1000, 1000000);
#else
  TFormula::AddVariable("pt");
  TFormula::AddVariable("eta");
  TFormula::AddVariable("phi");
  TFormula::AddVariable("energy");
  TFormula::AddVariable("d0");
  TFormula::AddVariable("dz");
  TFormula::AddVariable("ctgTheta");
#endif

  if(TFormula::Compile(buffer) != 0)
  {
    throw runtime_error("Invalid Long Formula.");
  }

  return 0;
}

//------------------------------------------------------------------------------

Double_t DelphesLongFormula::Eval(Double_t pt, Double_t eta, Double_t phi,
  Double_t energy, Double_t d0, Double_t dz, Double_t ctgTheta)
{
  Double_t x[7];

  x[0] = pt;
  x[1] = eta;
  x[2] = phi;
  x[3] = energy;
  x[4] = d0;
  x[5] = dz;
  x[6] = ctgTheta;

  return EvalPar(x);
}

//------------------------------------------------------------------------------
