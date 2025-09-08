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

 /** \class UnstablePropagator
  *
  *  Propagates charged unstable particles in magnetic field
  *  and updates coordinates of its daughters iteratively
  *
  *  \author M. Selvaggi - CERN
  *
  */

#include "modules/UnstablePropagator.h"

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

UnstablePropagator::UnstablePropagator() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

UnstablePropagator::~UnstablePropagator()
{
}

//------------------------------------------------------------------------------

void UnstablePropagator::Init()
{
  fRadius = GetDouble("Radius", 1.0);
  fRadius2 = fRadius * fRadius;
  fHalfLength = GetDouble("HalfLength", 3.0);
  fBz = GetDouble("Bz", 0.0);
  fLmin = GetDouble("Lmin", 1.0E-03);

  fDebug = false;
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

  fRadiusMax = GetDouble("RadiusMax", fRadius);
  fHalfLengthMax = GetDouble("HalfLengthMax", fHalfLength);
  // import array with output from filter/classifier module

  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();

}

//------------------------------------------------------------------------------

void UnstablePropagator::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void UnstablePropagator::Process()
{
  Candidate *candidate, *daughter;
  TLorentzVector particlePosition, particleMomentum;
  Double_t pt2, q;
  Double_t lof, x, y, z;

  fItInputArray->Reset();

  if (fDebug) cout<< "-------------   new event -----------------" << endl;

  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    particlePosition = candidate->Position;
    particleMomentum = candidate->Momentum;

    x = particlePosition.X() * 1.0E-3;
    y = particlePosition.Y() * 1.0E-3;
    z = particlePosition.Z() * 1.0E-3;
    pt2 = particleMomentum.Perp2();
    q = candidate->Charge;

    //if (fDebug) PrintPart("", candidate);

    // check that particle position is inside the cylinder
    if(TMath::Hypot(x, y) > fRadiusMax || TMath::Abs(z) > fHalfLengthMax)
    {
      continue;
    }

    if (TMath::Abs(q) < 1.0E-9 || TMath::Abs(fBz) < 1.0E-9)
    {
      continue;
    }

    if(pt2 < 1.0E-9)
    {
      continue;
    }

    // pass if particle already processed
    if (candidate->L > 1.0E-9)
    {
      continue;
    }

    std::vector < Int_t > daughters_indices = DaughterIndices(candidate);

    if (daughters_indices.size() == 0)
    {
      continue;
    }

    daughter = static_cast<Candidate *>(fInputArray->At(daughters_indices.at(0)));
    lof = FlightDistance(candidate, daughter) * 1.0E-3 ;

    //fLmin = 0.01;
    if (lof < fLmin)
    {
      continue;
    }

    if (fDebug) std::cout<<" -- lof: "<<lof<<", Lmin: "<<fLmin<<std::scientific<<std::endl;
    TString prefix = " -- ";
    ComputeChainFlightDistances(prefix, candidate);
    PropagateAndUpdateChain(prefix, candidate);

  }
}

//------------------------------------------------------------------------------

std::vector < Int_t > UnstablePropagator::DaughterIndices(Candidate *candidate)
{

  std::vector < Int_t > indices;

  Int_t d1 = candidate->D1;
  Int_t d2 = candidate->D2;

  Int_t maxd = max(d1, d2);
  Int_t mind = min(d1, d2);


  if (maxd < 0)
  {
    indices.clear();
  }
  else if (mind < 0)
  {
    indices.push_back(maxd);
  }
  else if (d1 > d2)
  {
    indices.push_back(d1);
    indices.push_back(d2);
  }
  else
  {
    for (Int_t i = d1; i <= d2; ++i)
    {
      indices.push_back(i);
    }
  }
  return indices;

}

//------------------------------------------------------------------------------

// returns flight distance in mm
Double_t UnstablePropagator::FlightDistance(Candidate *mother, Candidate *daughter)
{
  TVector3 vector = mother->Position.Vect() - daughter->Position.Vect();
  return vector.Mag() ;
}

//------------------------------------------------------------------------------

void UnstablePropagator::ComputeChainFlightDistances(TString prefix, Candidate *candidate)
{

  Candidate *daughter, *mother;
  mother = candidate;
  std::vector < Int_t > drange = DaughterIndices(mother);

  if (fDebug) cout<<prefix<<" computing chain flight distances"<<endl;
  if (fDebug) PrintPart(prefix, mother);
  prefix += " -- ";
  // check if particle already processed or if stable
  if (mother->L > 1.0E-9 || drange.size() == 0)
  //if (drange.size() == 0)
  {
    return;
  }
  else
  {
    daughter = static_cast<Candidate *>(fInputArray->At(drange.at(0)));
    mother->L = FlightDistance(mother, daughter);
    if (fDebug) cout<<prefix<<" flight distance: "<<mother->L<<endl;
    for(unsigned long i = 0; i < drange.size(); i++)
    {
      daughter = static_cast<Candidate *>(fInputArray->At(drange.at(i)));
      ComputeChainFlightDistances(prefix, daughter);
    }
  }
}

//------------------------------------------------------------------------------

void UnstablePropagator::PropagateAndUpdateChain(TString prefix, Candidate *candidate)
{
  Candidate *daughter, *mother;
  TLorentzVector updatedPosition;
  mother = candidate;
  std::vector < Int_t > drange = DaughterIndices(mother);

  //if (fDebug) cout<<prefix<<" propagating and updating chain, mother:"<<endl;
  if (fDebug) PrintPart(prefix, mother);
  //if (fDebug) cout<<mother->L<<","<<drange.size()<<endl;

  prefix += " --";

  // check if particle stable
  if (drange.size() == 0)
  {
    return;
  }
  else
  {
    updatedPosition = PropagatedPosition(mother);
    for(unsigned long i = 0; i < drange.size(); i++)
    {
      daughter = static_cast<Candidate *>(fInputArray->At(drange.at(i)));
      //  if (fDebug) cout<<prefix<<" propagating and updating chain, daughter:"<<endl;
      if (fDebug) PrintPart(prefix, daughter);
      daughter->Position = updatedPosition;
      //if (fDebug) cout<<prefix<<" propagated position: "<<daughter->Position.X()<<", "<<daughter->Position.Y()<<", "<<daughter->Position.Z()<<endl;
      PropagateAndUpdateChain(prefix, daughter);
    }
  }

}

//------------------------------------------------------------------------------

TLorentzVector UnstablePropagator::PropagatedPosition(Candidate *candidate)
{

  TLorentzVector particlePosition, particleMomentum, beamSpotPosition;
  Double_t px, py, pz, pt, e, q;
  Double_t x, y, z, t, r;
  Double_t x_c, y_c, phi_0;
  Double_t x_t, y_t, z_t, r_t, phi_t;
  Double_t gammam, omega;
  Double_t vz;
  Double_t tof, lof;

  const Double_t c_light = 2.99792458E8;

  particlePosition = candidate->Position;
  particleMomentum = candidate->Momentum;

  x = particlePosition.X() * 1.0E-3;
  y = particlePosition.Y() * 1.0E-3;
  z = particlePosition.Z() * 1.0E-3;

  q = candidate->Charge;

  px = particleMomentum.Px();
  py = particleMomentum.Py();
  pz = particleMomentum.Pz();
  pt = particleMomentum.Pt();
  e = particleMomentum.E();

  // propagation flight and time of flight
  lof = candidate->L * 1.0E-3; // in meters
  tof = lof/(particleMomentum.Beta() * c_light); // in seconds

  //if (fDebug) cout << "propagating from : "<<x<<", "<<y<<", "<<z<<",  lof:"<<lof<<", tof: "<<tof<<endl;

  if(TMath::Hypot(x, y) > fRadius || TMath::Abs(z) > fHalfLength)
  {
    return particlePosition;
  }

  // neutral propagation
  else if(TMath::Abs(q) < 1.0E-9 || TMath::Abs(fBz) < 1.0E-9)
  {
    // solve pt2*t^2 + 2*(px*x + py*y)*t - (fRadius2 - x*x - y*y) = 0

    /*
    tmp = px * y - py * x;
    t_r = (TMath::Sqrt(pt2 * fRadius2 - tmp * tmp) - px * x - py * y) / pt2;

    t_z = (TMath::Sign(fHalfLength, pz) - z) / pz;

    t = TMath::Min(t_r, t_z);
    */

    // TODO: check that l and t within cilinder
    // TODO: convert time properly in proper units
    t =  c_light * tof / e ;
    x_t = x + px * t;
    y_t = y + py * t;
    z_t = z + pz * t;
    // Double_t l = TMath::Sqrt((x_t - x) * (x_t - x) + (y_t - y) * (y_t - y) + (z_t - z) * (z_t - z));
    //if (fDebug) cout << "propagated neutral length: "<<l<<endl;
    particlePosition.SetXYZT(x_t * 1.0E3, y_t * 1.0E3, z_t * 1.0E3, particlePosition.T() + t * e * 1.0E3);
  }
  else
  {

    // 1. initial transverse momentum p_{T0}: Part->pt
    //    initial transverse momentum direction phi_0 = -atan(p_{X0} / p_{Y0})
    //    relativistic gamma: gamma = E / mc^2; gammam = gamma * m
    //    gyration frequency omega = q * Bz / (gammam)
    //    helix radius r = p_{T0} / (omega * gammam)

    gammam = e * 1.0E9 / (c_light * c_light); // gammam in [eV/c^2]
    omega = q * fBz / gammam; // omega is here in [89875518/s]
    r = pt / (q * fBz) * 1.0E9 / c_light; // in [m]

    phi_0 = TMath::ATan2(py, px); // [rad] in [-pi, pi]

    // 2. helix axis coordinates
    x_c = x + r * TMath::Sin(phi_0);
    y_c = y - r * TMath::Cos(phi_0);

    vz = pz * c_light / e;

    // 3. position in terms of x(t), y(t), z(t)

    t = tof; // in seconds
    phi_t = phi_0 - omega * t;
    x_t = x_c - r * TMath::Sin(phi_t);
    y_t = y_c + r * TMath::Cos(phi_t);
    z_t = z + vz * t;
    r_t = TMath::Hypot(x_t, y_t);

    // lenght of the path from production to tracker
    // l = t * TMath::Hypot(vz, r * omega);
    //if (fDebug) cout << "propagated to: X: "<<x_t<<" Y: "<<y_t<<" Z: "<<z_t<<"  length: "<<l<<endl;

    if(r_t > 0.0)
    {
      particlePosition.SetXYZT(x_t * 1.0E3, y_t * 1.0E3, z_t * 1.0E3, particlePosition.T() + t * c_light * 1.0E3);
    }
  }

  return particlePosition;
}

//------------------------------------------------------------------------------

void UnstablePropagator::PrintPart(TString prefix, Candidate *candidate)
{
  std::cout.precision(6);
  std::cout<<prefix;
  std::cout<<Index(candidate)<<", PID:"<<candidate->PID <<", Q:"<<std::scientific;
  std::cout<<candidate->Charge <<", Status:"<<candidate->Status<<", E:"<<std::scientific;
  std::cout<<candidate->Momentum.E() <<", Eta:"<<candidate->Momentum.Eta()<<", Phi:"<<std::scientific;
  std::cout<<candidate->Momentum.Phi() <<", X:"<<candidate->Position.X()<<", Y:"<<std::scientific;
  std::cout<<candidate->Position.Y() <<", Z:"<<candidate->Position.Z()<<", D1:"<<std::scientific;
  std::cout<<candidate->D1 <<", D2:"<<candidate->D2;
  std::cout<<", L:"<<candidate->L<<std::scientific<<std::endl;
}

//------------------------------------------------------------------------------

Int_t UnstablePropagator::Index(Candidate *particle)
{
  /*Int_t i=-1;
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    i++;
    if(candidate->GetUniqueID() == particle->GetUniqueID())
    {
      break;
    }
  }
  */
  Int_t j=-1;
  for(Int_t i = 0; i < fInputArray->GetEntriesFast(); i++)
  {
    j = i;
    if(fInputArray->At(i)->GetUniqueID() == particle->GetUniqueID())
    {
      break;
    }
  }
  return j;
}
