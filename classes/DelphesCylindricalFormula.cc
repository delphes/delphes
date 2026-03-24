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

#include <stdexcept>

//------------------------------------------------------------------------------

DelphesCylindricalFormula::DelphesCylindricalFormula() :
  TFormula()
{
}

//------------------------------------------------------------------------------

DelphesCylindricalFormula::DelphesCylindricalFormula(std::string_view /*name*/, std::string_view /*expression*/) :
  TFormula()
{
}

//------------------------------------------------------------------------------

Int_t DelphesCylindricalFormula::Compile(std::string_view expression)
{
  TString buffer;
  const char *it;
  for(it = expression.data(); *it; ++it)
  {
    if(*it == ' ' || *it == '\t' || *it == '\r' || *it == '\n' || *it == '\\') continue;
    buffer.Append(*it);
  }
  buffer.ReplaceAll("r", "x");
  buffer.ReplaceAll("phi", "y");
  buffer.ReplaceAll("z", "z");
  if(TFormula::Compile(buffer) != 0)
    throw std::runtime_error("Invalid formula.");
  return 0;
}

//------------------------------------------------------------------------------

double DelphesCylindricalFormula::Eval(double r, double phi, double z) const
{
  double x[3] = {r, phi, z};
  return EvalPar(x);
}

//------------------------------------------------------------------------------
