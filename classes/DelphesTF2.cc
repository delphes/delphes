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
