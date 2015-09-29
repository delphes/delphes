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


/** \class PhotonConversions
 *
 *
 * Converts photons into e+ e- pairs according to mass ditribution in the detector.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/PhotonConversions.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesCylindricalFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TF1.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

PhotonConversions::PhotonConversions() :
  fItInputArray(0), fConversionMap(0), fDecayXsec(0)
{
  fDecayXsec = new TF1;
  fConversionMap = new DelphesCylindricalFormula;
}

//------------------------------------------------------------------------------

PhotonConversions::~PhotonConversions()
{
}

//------------------------------------------------------------------------------

void PhotonConversions::Init()
{
  fRadius = GetDouble("Radius", 1.0);
  fRadius2 = fRadius*fRadius;

  fHalfLength = GetDouble("HalfLength", 3.0);

  fEtaMax = GetDouble("EtaMax", 5.0);
  fEtaMin = GetDouble("EtaMin", 2.0);

  fStep = GetDouble("Step", 0.1); // in meters

  fConversionMap->Compile(GetString("ConversionMap", "0.0"));

#if  ROOT_VERSION_CODE < ROOT_VERSION(6,04,00)
  fDecayXsec->Compile("1.0 - 4.0/3.0 * x * (1.0 - x)");
#else
  fDecayXsec->GetFormula()->Compile("1.0 - 4.0/3.0 * x * (1.0 - x)");
#endif
  fDecayXsec->SetRange(0.0, 1.0);

  // import array with output from filter/classifier module

  fInputArray = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output arrays

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void PhotonConversions::Finish()
{
  if(fItInputArray) delete fItInputArray;
  if(fDecayXsec) delete fDecayXsec;
  if(fConversionMap) delete fConversionMap;
}

//------------------------------------------------------------------------------

void PhotonConversions::Process()
{
  Candidate *candidate, *ep, *em;
  TLorentzVector candidatePosition, candidateMomentum;
  TVector3 pos_i;
  Double_t px, py, pz, pt, pt2, e, eta, phi;
  Double_t x, y, z, t;
  Double_t x_t, y_t, z_t, r_t;
  Double_t x_i, y_i, z_i, r_i, phi_i;
  Double_t dt, t1, t2, t3, t4;
  Double_t tmp, discr, discr2;
  Int_t nsteps, i;
  Double_t rate, p_conv, x1, x2;
  Bool_t converted;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {

    if(candidate->PID != 22)
    {
      fOutputArray->Add(candidate);
    }
    else
    {
      candidatePosition = candidate->Position;
      candidateMomentum = candidate->Momentum;
      x = candidatePosition.X()*1.0E-3;
      y = candidatePosition.Y()*1.0E-3;
      z = candidatePosition.Z()*1.0E-3;

      // check that particle position is inside the cylinder
      if(TMath::Hypot(x, y) > fRadius || TMath::Abs(z) > fHalfLength) continue;

      px = candidateMomentum.Px();
      py = candidateMomentum.Py();
      pz = candidateMomentum.Pz();
      pt = candidateMomentum.Pt();
      pt2 = candidateMomentum.Perp2();
      eta = candidateMomentum.Eta();
      phi = candidateMomentum.Phi();
      e = candidateMomentum.E();

      if(eta < fEtaMin || eta > fEtaMax) continue;

      // solve pt2*t^2 + 2*(px*x + py*y)*t - (fRadius2 - x*x - y*y) = 0
      tmp = px*y - py*x;
      discr2 = pt2*fRadius2 - tmp*tmp;

      if(discr2 < 0.0)
      {
        // no solutions
        continue;
      }

      tmp = px*x + py*y;
      discr = TMath::Sqrt(discr2);
      t1 = (-tmp + discr)/pt2;
      t2 = (-tmp - discr)/pt2;
      t = (t1 < 0.0) ? t2 : t1;

      z_t = z + pz*t;
      if(TMath::Abs(z_t) > fHalfLength)
      {
        t3 = (+fHalfLength - z) / pz;
        t4 = (-fHalfLength - z) / pz;
        t = (t3 < 0.0) ? t4 : t3;
      }

      // final position
      x_t = x + px*t;
      y_t = y + py*t;
      z_t = z + pz*t;

      r_t = TMath::Sqrt(x_t*x_t + y_t*y_t + z_t*z_t);


      // here starts conversion code
      nsteps = Int_t(r_t/fStep);

      x_i = x;
      y_i = y;
      z_i = z;

      dt = t/nsteps;

      converted = false;

      for(i = 0; i < nsteps; ++i)
      {
        x_i += px*dt;
        y_i += py*dt;
        z_i += pz*dt;
        pos_i.SetXYZ(x_i,y_i,z_i);

        // convert photon position into cylindrical coordinates, cylindrical r,phi,z !!

        r_i = TMath::Sqrt(x_i*x_i + y_i*y_i);
        phi_i = pos_i.Phi();

        // read conversion rate/meter from card
        rate = fConversionMap->Eval(r_i, phi_i, z_i);

        // convert into conversion probability
        p_conv = 1 - TMath::Exp(-7.0/9.0*fStep*rate);

        // case conversion occurs
        if(gRandom->Uniform() < p_conv)
        {
          converted = true;

          // generate x1 and x2, the fraction of the photon energy taken resp. by e+ and e-
          x1 = fDecayXsec->GetRandom();
          x2 = 1 - x1;

          ep = static_cast<Candidate*>(candidate->Clone());
          em = static_cast<Candidate*>(candidate->Clone());

          ep->Position.SetXYZT(x_i*1.0E3, y_i*1.0E3, z_i*1.0E3, candidatePosition.T() + nsteps*dt*e*1.0E3);
          em->Position.SetXYZT(x_i*1.0E3, y_i*1.0E3, z_i*1.0E3, candidatePosition.T() + nsteps*dt*e*1.0E3);

          ep->Momentum.SetPtEtaPhiE(x1*pt, eta, phi, x1*e);
          em->Momentum.SetPtEtaPhiE(x2*pt, eta, phi, x2*e);

          ep->PID = -11;
          em->PID = 11;

          ep->Charge = 1.0;
          em->Charge = -1.0;

          ep->IsFromConversion = 1;
          em->IsFromConversion = 1;

          fOutputArray->Add(em);
          fOutputArray->Add(ep);

          break;
        }
      }
      if(!converted) fOutputArray->Add(candidate);
    }
  }
}

//------------------------------------------------------------------------------

