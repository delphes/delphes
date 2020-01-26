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

#ifndef EnergyLoss_h
#define EnergyLoss_h

/** \class EnergyLoss
 *
 *  This module computes the charged energy loss according to the active material properties.
 *  The energy loss is simulated with a Landau convoluted by a Gaussian.
 *
 *  \author M. Selvaggi - CERN
 *
 */

#include "classes/DelphesModule.h"
#include <map>
#include "TF1.h"

class TIterator;
class TObjArray;
class DelphesFormula;

class EnergyLoss: public DelphesModule
{
public:
  EnergyLoss();
  ~EnergyLoss();

  void Init();
  void Process();
  void Finish();

private:

  Double_t   fActiveFraction;
  Double_t   fThickness;

  Double_t   fResolution;
  Double_t   fTruncatedMeanFraction;

  // material parameters
  Double_t   fZ;
  Double_t   fA;
  Double_t   fRho;

  Double_t   fAa;
  Double_t   fM;
  Double_t   fX0;
  Double_t   fX1;
  Double_t   fI;
  Double_t   fC0;

  // this function computes corrections due to polarisation of the material
  Double_t Deltaf(Double_t c0, Double_t a, Double_t m, Double_t x0, Double_t x1, Double_t beta, Double_t gamma);

  // this function computes truncated of dEdx measurements
  Double_t TruncatedMean(std::vector<Double_t> elosses, Double_t truncFrac);

  std::vector<TIterator *> fInputList; //!

  ClassDef(EnergyLoss, 1)

};

#endif
