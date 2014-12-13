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

