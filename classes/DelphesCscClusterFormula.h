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

#ifndef DelphesCscClusterFormula_h
#define DelphesCscClusterFormula_h

#include "TFormula.h"

class Candidate;

class DelphesCscClusterFormula: public TFormula
{
public:
  DelphesCscClusterFormula();

  DelphesCscClusterFormula(const char *name, const char *expression);

  ~DelphesCscClusterFormula();

  Int_t Compile(const char *expression);

  Double_t Eval(Double_t deltaR, Double_t deltaZ = 0, Double_t Ehad = 0, Double_t Eem = 0);
};

#endif /* DelphesCscClusterFormula_h */
