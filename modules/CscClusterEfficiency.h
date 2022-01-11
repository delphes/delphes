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

#ifndef CscClusterEfficiency_h
#define CscClusterEfficiency_h



/** \class CscClusterEfficiency
 *
 *  This module is specific to the CMS paper searching for neutral LLPs in the CMS endcap muon detectors: https://arxiv.org/abs/2107.04838
 *  It is implemented based on the ClusterEfficiency parameterization function provided in the HEPData entry of the paper: https://www.hepdata.net/record/104408
 *
 *  \author Christina Wang
 *
 */
#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;
class DelphesCscClusterFormula;

class CscClusterEfficiency: public DelphesModule
{
public:
  CscClusterEfficiency();
  ~CscClusterEfficiency();

  void Init();
  void Process();
  void Finish();

private:
  DelphesCscClusterFormula *fFormula; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(CscClusterEfficiency, 1)
};

#endif
