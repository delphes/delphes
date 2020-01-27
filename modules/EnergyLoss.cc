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

 /** \class EnergyLoss
  *
  *  This module computes the charged energy loss according to the active material properties.
  *  The energy loss is simulated with a Landau convoluted by a Gaussian.
  *
  *  \author M. Selvaggi - CERN
  *
  */
#include "modules/EnergyLoss.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

EnergyLoss::EnergyLoss()
{
}

//------------------------------------------------------------------------------

EnergyLoss::~EnergyLoss()
{
}

//------------------------------------------------------------------------------

void EnergyLoss::Init()
{

  fActiveFraction              = GetDouble("ActiveFraction", 0.002); // active fraction of the detector
  fThickness                   = GetDouble("Thickness", 200E-6); // active detector thickness
  fResolution                  = GetDouble("Resolution", 0.4); // 0 - perfect Landau energy loss (0.15 gives good agreement with CMS pixel detector)
  fTruncatedMeanFraction       = GetDouble("TruncatedMeanFraction", 0.5); // fraction of measurements to ignore when computing mean

  // active material properties (cf. http://pdg.lbl.gov/2014/AtomicNuclearProperties/properties8.dat)
  fZ    =  GetDouble("Z", 14.);
  fA    =  GetDouble("A",  28.0855); // in g/mol
  fRho  =  GetDouble("rho", 2.329); // in g/cm3
  fAa   =  GetDouble("a", 0.1492);
  fM    =  GetDouble("m", 3.2546);
  fX0   =  GetDouble("x0", 0.2015);
  fX1   =  GetDouble("x1", 2.8716);
  fI    =  GetDouble("I", 173.0); // mean excitation potential in (eV)
  fC0   =  GetDouble("c0", 4.4355);

  // import arrays with output from other modules

  ExRootConfParam param = GetParam("InputArray");
  Long_t i, size;
  const TObjArray *array;
  TIterator *iterator;

  size = param.GetSize();
  for(i = 0; i < size; ++i)
  {
    array = ImportArray(param[i].GetString());
    iterator = array->MakeIterator();

    fInputList.push_back(iterator);
  }

}

//------------------------------------------------------------------------------

void EnergyLoss::Finish()
{
  vector<TIterator *>::iterator itInputList;
  TIterator *iterator;

  for(itInputList = fInputList.begin(); itInputList != fInputList.end(); ++itInputList)
  {
    iterator = *itInputList;
    if(iterator) delete iterator;
  }

}

//------------------------------------------------------------------------------

void EnergyLoss::Process()
{
  Candidate *candidate;
  vector<TIterator *>::iterator itInputList;
  TIterator *iterator;

  Double_t beta, gamma, charge;
  Double_t kappa, chi, me, I, Wmax, delta, avdE, dP, dx, L, dE, dEdx, res;
  Double_t eloss_truncmean;

  Int_t nhits;
  vector<Double_t> elosses;

  //cout<<"---------------- new event -------------------"<<endl;


  // loop over all input arrays
  for(itInputList = fInputList.begin(); itInputList != fInputList.end(); ++itInputList)
  {
    iterator = *itInputList;

    // loop over all candidates
    iterator->Reset();
    while((candidate = static_cast<Candidate *>(iterator->Next())))
    {
      //cout<<"    ---------------- new candidate -------------------"<<endl;
      const TLorentzVector &candidateMomentum = candidate->Momentum;

      beta      = candidateMomentum.Beta();
      gamma     = candidateMomentum.Gamma();
      charge    = TMath::Abs(candidate->Charge);


      // length of the track normalized by the fraction of active material and the charge collection efficiency in the tracker (in cm)
      //dx = candidate->L * fActiveFraction * 0.1;
      // amount of material in one sensor (converted in cm)
      dx = fThickness * 100.;
      // path length in cm
      L  = candidate->L * 0.1;


     // compute number of hits as path length over active length
      nhits = Int_t(L*fActiveFraction/dx);


      //beta = 0.999945;
      //gamma = 95.6446;
      //charge = 1.;
      //nhits = 100;      

      //cout<<L<<","<<fActiveFraction<<","<<dx<<","<<nhits<<endl;

      kappa = 2*0.1535*TMath::Abs(charge)*TMath::Abs(charge)*fZ*fRho*dx/(fA*beta*beta); //energy loss in MeV

      chi = 0.5*kappa;
      me = 0.510998; // electron mass in MeV, need
      I = fI*1e-6; // convert I in MeV

      // fixme: max energy transfer wrong for electrons
      Wmax = 2*me*beta*beta*gamma*gamma; // this is not valid for electrons

      delta = Deltaf(fC0, fAa, fM, fX0, fX1, beta, gamma);

      // Bethe-Bloch energy loss in MeV (not used here)
      avdE  = kappa*( TMath::Log(Wmax/I) - beta*beta - delta/2);

      // most probable energy (MPV) loss for Landau in a single layer
      dP = chi*( TMath::Log(Wmax/I) + TMath::Log(chi/I) + 0.2 - beta*beta - delta);

      //cout<<"L: "<<L<<", PT: "<<candidateMomentum.Pt()<<", Eta: "<<candidateMomentum.Eta()<<", Phi: "<< candidateMomentum.Phi()<<endl;
      //cout<<"Nhits: "<<nhits<<", dx: "<<dx<<", Charge: "<<charge<<", Beta: "<< beta<<", Gamma: "<<gamma<<", PT: "<<candidateMomentum.Pt()<<endl;
      //cout<<x<<","<<kappa<<endl;
      //cout<<"    Wmax: "<<Wmax<<", Chi: "<<chi<<", delta: "<<delta<<", DeDx: "<<avdE<<", DeltaP: "<<dP<<endl;

      // simulate Nhits energy loss measurements
      elosses.clear();
      for (Int_t j=0; j<nhits; j++){
        // compute total energy loss in MeV predicted by a Landau
        dE = gRandom->Landau(dP,chi); // this is the total energy loss in MeV predicted by a Landau

        // convert resolution given in Mev/cm into absolute for this sensor
        res = fResolution*dx;

        // apply additionnal gaussian smearing
        dE = gRandom->Gaus(dE,res);
        elosses.push_back(dE);
      }

      sort (elosses.begin(), elosses.end());
      eloss_truncmean = TruncatedMean(elosses, fTruncatedMeanFraction);


      dEdx = dx > 0 ? eloss_truncmean/dx : -1. ;

      // store computed dEdx in MeV/cm
      candidate->DeDx = dEdx;

      // add dedx also in Muons in electrons classes in treeWriter
      // fix electrons here
      // think whether any relevance for hits


      //cout<<"    eloss: "<<dE<<", dx: "<<dx<<", dEdx: "<<dEdx<<endl;
    }
  }

}

//------------------------------------------------------------------------------

// formula Taken from Leo (2.30) pg. 26
Double_t EnergyLoss::Deltaf(Double_t c0, Double_t a, Double_t m, Double_t x0, Double_t x1, Double_t beta, Double_t gamma)
{
   Double_t x= TMath::Log10(beta*gamma);
   Double_t delta = 0.;

   //cout<<x<<","<<x0<<","<<x1<<","<<endl;

   if (x < x0)
       delta = 0.;
   if (x >= x0 && x< x1)
       delta = 4.6052*x - c0 + a*TMath::Power(x1 - x,m);
   if  (x> x1)
       delta = 4.6052*x - c0;

   return delta;
}

//------------------------------------------------------------------------------
Double_t EnergyLoss::TruncatedMean(std::vector<Double_t> elosses, Double_t truncFrac)
{
     Int_t new_size = Int_t( elosses.size() * (1 - truncFrac));

     // remove outliers and re-compute mean
     elosses.resize(new_size);
     return accumulate( elosses.begin(), elosses.end(), 0.0)/elosses.size();
}
