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
 *  Subtract pile-up contribution from tracks.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
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
  Double_t fActiveFraction;
  Double_t fChargeCollectionEfficiency;
  Double_t fResolution;  

  // material parameters
  Double_t   fZ;   
  Double_t   fA;   
  Double_t   frho; 

  Double_t   fa;   
  Double_t   fm;   
  Double_t   fx0;  
  Double_t   fx1;  
  Double_t   fI;  
  Double_t   fc0;  

  // this function computes corrections due to polarisation of the material
  Double_t Deltaf(Double_t c0, Double_t a, Double_t m, Double_t x0, Double_t x1, Double_t beta, Double_t gamma);

  std::vector<TIterator *> fInputList; //!

  ClassDef(EnergyLoss, 1)

};

#endif
