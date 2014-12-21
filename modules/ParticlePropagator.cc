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


/** \class ParticlePropagator
 *
 *  Propagates charged and neutral particles
 *  from a given vertex to a cylinder defined by its radius,
 *  its half-length, centered at (0,0,0) and with its axis
 *  oriented along the z-axis.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/ParticlePropagator.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

ParticlePropagator::ParticlePropagator() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

ParticlePropagator::~ParticlePropagator()
{
}

//------------------------------------------------------------------------------

void ParticlePropagator::Init()
{
  fRadius = GetDouble("Radius", 1.0);
  fRadius2 = fRadius*fRadius;
  fHalfLength = GetDouble("HalfLength", 3.0);
  fBz = GetDouble("Bz", 0.0);
  if(fRadius < 1.0E-2)
  {
    cout << "ERROR: magnetic field radius is too low\n";
    return;
  }
  if(fHalfLength < 1.0E-2)
  {
    cout << "ERROR: magnetic field length is too low\n";
    return;
  }

  // import array with output from filter/classifier module

  fInputArray = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output arrays

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
  fChargedHadronOutputArray = ExportArray(GetString("ChargedHadronOutputArray", "chargedHadrons"));
  fElectronOutputArray = ExportArray(GetString("ElectronOutputArray", "electrons"));
  fMuonOutputArray = ExportArray(GetString("MuonOutputArray", "muons"));
}

//------------------------------------------------------------------------------

void ParticlePropagator::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void ParticlePropagator::Process()
{
  Candidate *candidate, *mother;
  TLorentzVector candidatePosition, candidateMomentum;
  Double_t px, py, pz, pt, pt2, e, q;
  Double_t x, y, z, t, r, phi;
  Double_t x_c, y_c, r_c, phi_c, phi_0;
  Double_t x_t, y_t, z_t, r_t;
  Double_t t1, t2, t3, t4, t5, t6;
  Double_t t_z, t_r, t_ra, t_rb;
  Double_t tmp, discr, discr2;
  Double_t delta, gammam, omega, asinrho;
  Double_t rcu, rc2, dxy, xd, yd, zd;

  const Double_t c_light = 2.99792458E8;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    candidatePosition = candidate->Position;
    candidateMomentum = candidate->Momentum;
    x = candidatePosition.X()*1.0E-3;
    y = candidatePosition.Y()*1.0E-3;
    z = candidatePosition.Z()*1.0E-3;
    q = candidate->Charge;

    // check that particle position is inside the cylinder
    if(TMath::Hypot(x, y) > fRadius || TMath::Abs(z) > fHalfLength)
    {
      continue;
    }

    px = candidateMomentum.Px();
    py = candidateMomentum.Py();
    pz = candidateMomentum.Pz();
    pt = candidateMomentum.Pt();
    pt2 = candidateMomentum.Perp2();
    e = candidateMomentum.E();

    if(pt2 < 1.0E-9)
    {
      continue;
    }

    if(TMath::Abs(q) < 1.0E-9 || TMath::Abs(fBz) < 1.0E-9)
    {
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

      x_t = x + px*t;
      y_t = y + py*t;
      z_t = z + pz*t;

      mother = candidate;
      candidate = static_cast<Candidate*>(candidate->Clone());

      candidate->Position.SetXYZT(x_t*1.0E3, y_t*1.0E3, z_t*1.0E3, candidatePosition.T() + t*e*1.0E3);

      candidate->Momentum = candidateMomentum;
      candidate->AddCandidate(mother);

      fOutputArray->Add(candidate);
      if(TMath::Abs(q) > 1.0E-9)
      {
        switch(TMath::Abs(candidate->PID))
        {
          case 11:
            fElectronOutputArray->Add(candidate);
            break;
          case 13:
            fMuonOutputArray->Add(candidate);
            break;
          default:
            fChargedHadronOutputArray->Add(candidate);
        }
      }
    }
    else
    {

      // 1.  initial transverse momentum p_{T0}: Part->pt
      //     initial transverse momentum direction phi_0 = -atan(p_X0/p_Y0)
      //     relativistic gamma: gamma = E/mc^2; gammam = gamma * m
      //     gyration frequency omega = q/(gamma m) fBz
      //     helix radius r = p_{T0} / (omega gamma m)

      gammam = e*1.0E9 / (c_light*c_light);      // gammam in [eV/c^2]
      omega = q * fBz / (gammam);                // omega is here in [89875518/s]
      r = pt / (q * fBz) * 1.0E9/c_light;        // in [m]

      phi_0 = TMath::ATan2(py, px); // [rad] in [-pi, pi]

      // 2. helix axis coordinates
      x_c = x + r*TMath::Sin(phi_0);
      y_c = y - r*TMath::Cos(phi_0);
      r_c = TMath::Hypot(x_c, y_c);
      phi_c = TMath::ATan2(y_c, x_c);
      phi = phi_c;
      if(x_c < 0.0) phi += TMath::Pi();

      rcu = TMath::Abs(r);
      rc2 = r_c*r_c;

      // calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
      xd = x_c*x_c*x_c - x_c*rcu*r_c + x_c*y_c*y_c;
      xd = (rc2 > 0.0) ? xd / rc2 : -999;
      yd = y_c*(-rcu*r_c + rc2);
      yd = (rc2 > 0.0) ? yd / rc2 : -999;
      zd = z + (TMath::Sqrt(xd*xd + yd*yd) - TMath::Sqrt(x*x + y*y))*pz/pt;

      // calculate impact paramater
      dxy = (xd*py - yd*px)/pt;

      // 3. time evaluation t = TMath::Min(t_r, t_z)
      //    t_r : time to exit from the sides
      //    t_z : time to exit from the front or the back
      t_r = 0.0; // in [ns]
      int sign_pz = (pz > 0.0) ? 1 : -1;
      if(pz == 0.0) t_z = 1.0E99;
      else t_z = gammam / (pz*1.0E9/c_light) * (-z + fHalfLength*sign_pz);

      if(r_c + TMath::Abs(r)  < fRadius)
      {
        // helix does not cross the cylinder sides
        t = t_z;
      }
      else
      {
        asinrho = TMath::ASin( (fRadius*fRadius - r_c*r_c - r*r) / (2*TMath::Abs(r)*r_c)  );
        delta = phi_0 - phi;
        if(delta <-TMath::Pi()) delta += 2*TMath::Pi();
        if(delta > TMath::Pi()) delta -= 2*TMath::Pi();
        t1 = (delta + asinrho) / omega;
        t2 = (delta + TMath::Pi() - asinrho) / omega;
        t3 = (delta + TMath::Pi() + asinrho) / omega;
        t4 = (delta - asinrho) / omega;
        t5 = (delta - TMath::Pi() - asinrho) / omega;
        t6 = (delta - TMath::Pi() + asinrho) / omega;

        if(t1 < 0.0) t1 = 1.0E99;
        if(t2 < 0.0) t2 = 1.0E99;
        if(t3 < 0.0) t3 = 1.0E99;
        if(t4 < 0.0) t4 = 1.0E99;
        if(t5 < 0.0) t5 = 1.0E99;
        if(t6 < 0.0) t6 = 1.0E99;

        t_ra = TMath::Min(t1, TMath::Min(t2, t3));
        t_rb = TMath::Min(t4, TMath::Min(t5, t6));
        t_r = TMath::Min(t_ra, t_rb);
        t = TMath::Min(t_r, t_z);
      }

      // 4. position in terms of x(t), y(t), z(t)
      x_t = x_c + r * TMath::Sin(omega * t - phi_0);
      y_t = y_c + r * TMath::Cos(omega * t - phi_0);
      z_t = z + pz*1.0E9 / c_light / gammam * t;
      r_t = TMath::Hypot(x_t, y_t);

      if(r_t > 0.0)
      {
        mother = candidate;
        candidate = static_cast<Candidate*>(candidate->Clone());

        candidate->Position.SetXYZT(x_t*1.0E3, y_t*1.0E3, z_t*1.0E3, candidatePosition.T() + t*c_light*1.0E3);

        candidate->Momentum = candidateMomentum;
        candidate->Dxy = dxy*1.0E3;
        candidate->Xd = xd*1.0E3;
        candidate->Yd = yd*1.0E3;
        candidate->Zd = zd*1.0E3;

        candidate->AddCandidate(mother);

        fOutputArray->Add(candidate);
        switch(TMath::Abs(candidate->PID))
        {
          case 11:
            fElectronOutputArray->Add(candidate);
            break;
          case 13:
            fMuonOutputArray->Add(candidate);
            break;
          default:
            fChargedHadronOutputArray->Add(candidate);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

