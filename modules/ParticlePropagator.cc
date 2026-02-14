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
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

void ParticlePropagator::Init()
{
  fRadius = GetDouble("Radius", 1.0);
  fRadius2 = fRadius * fRadius;
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

  fRadiusMax = GetDouble("RadiusMax", fRadius);
  fHalfLengthMax = GetDouble("HalfLengthMax", fHalfLength);

  // import array with output from filter/classifier module
  ImportArray(GetString("InputArray", "Delphes/stableParticles"), fInputArray);
  try
  { // import beamspot
    ImportArray(GetString("BeamSpotInputArray", "BeamSpotFilter/beamSpotParticle"), fBeamSpotInputArray);
  }
  catch(runtime_error &e)
  {
  }

  // create output arrays
  ExportArray(fOutputArray, GetString("OutputArray", "stableParticles"));
  ExportArray(fNeutralOutputArray, GetString("NeutralOutputArray", "neutralParticles"));
  ExportArray(fChargedHadronOutputArray, GetString("ChargedHadronOutputArray", "chargedHadrons"));
  ExportArray(fElectronOutputArray, GetString("ElectronOutputArray", "electrons"));
  ExportArray(fMuonOutputArray, GetString("MuonOutputArray", "muons"));
}

//------------------------------------------------------------------------------

void ParticlePropagator::Finish()
{
}

//------------------------------------------------------------------------------

void ParticlePropagator::Process()
{
  Double_t px, py, pz, pt, pt2, e, q;
  Double_t x, y, z, t, r;
  Double_t x_c, y_c, r_c, phi_0;
  Double_t x_t, y_t, z_t, r_t, phi_t;
  Double_t t_r, t_z;
  Double_t tmp;
  Double_t gammam, omega;
  Double_t xd, yd, zd;
  Double_t l, d0, dz, ctgTheta, alpha;
  Double_t bsx, bsy, bsz;
  Double_t td, pio, phid, vz;

  const Double_t c_light = 2.99792458E8;

  fOutputArray->clear();
  fNeutralOutputArray->clear();
  fChargedHadronOutputArray->clear();
  fElectronOutputArray->clear();
  fMuonOutputArray->clear();

  ROOT::Math::XYZTVector beamSpotPosition;
  if(fBeamSpotInputArray && !fBeamSpotInputArray->empty())
  {
    const auto &beamSpotCandidate = fBeamSpotInputArray->at(0);
    beamSpotPosition = beamSpotCandidate.Position;
  }

  for(auto &candidate : *fInputArray) //TODO: ensure a const-qualified version cannot be used
  {
    Candidate *particle = nullptr;
    if(candidate.GetCandidates().empty())
    {
      particle = &candidate;
    }
    else
    {
      particle = static_cast<Candidate *>(candidate.GetCandidates().at(0));
    }

    const auto &particlePosition = particle->Position;
    const auto &particleMomentum = particle->Momentum;

    x = particlePosition.X() * 1.0E-3;
    y = particlePosition.Y() * 1.0E-3;
    z = particlePosition.Z() * 1.0E-3;

    bsx = beamSpotPosition.X() * 1.0E-3;
    bsy = beamSpotPosition.Y() * 1.0E-3;
    bsz = beamSpotPosition.Z() * 1.0E-3;

    q = particle->Charge;

    // check that particle position is inside the cylinder
    if(TMath::Hypot(x, y) > fRadiusMax || TMath::Abs(z) > fHalfLengthMax)
    {
      continue;
    }

    px = particleMomentum.Px();
    py = particleMomentum.Py();
    pz = particleMomentum.Pz();
    pt = particleMomentum.Pt();
    pt2 = particleMomentum.Perp2();
    e = particleMomentum.E();

    if(pt2 < 1.0E-9)
    {
      continue;
    }

    if(TMath::Hypot(x, y) > fRadius || TMath::Abs(z) > fHalfLength)
    {
      auto new_candidate = candidate;

      new_candidate.InitialPosition = particlePosition;
      new_candidate.Position = particlePosition;
      new_candidate.L = 0.0;

      new_candidate.Momentum = particleMomentum;
      new_candidate.AddCandidate(&candidate); // preserve parentage

      fOutputArray->emplace_back(new_candidate);
    }
    else if(TMath::Abs(q) < 1.0E-9 || TMath::Abs(fBz) < 1.0E-9)
    {
      // solve pt2*t^2 + 2*(px*x + py*y)*t - (fRadius2 - x*x - y*y) = 0
      tmp = px * y - py * x;
      t_r = (TMath::Sqrt(pt2 * fRadius2 - tmp * tmp) - px * x - py * y) / pt2;

      t_z = (TMath::Sign(fHalfLength, pz) - z) / pz;

      t = TMath::Min(t_r, t_z);

      x_t = x + px * t;
      y_t = y + py * t;
      z_t = z + pz * t;

      l = TMath::Sqrt((x_t - x) * (x_t - x) + (y_t - y) * (y_t - y) + (z_t - z) * (z_t - z));

      auto new_candidate = candidate;

      new_candidate.InitialPosition = particlePosition;
      new_candidate.Position.SetXYZT(x_t * 1.0E3, y_t * 1.0E3, z_t * 1.0E3, particlePosition.T() + t * e * 1.0E3);
      new_candidate.L = l * 1.0E3;

      new_candidate.Momentum = particleMomentum;
      new_candidate.AddCandidate(&candidate); // preserve parentage

      fOutputArray->emplace_back(new_candidate);

      if(TMath::Abs(q) > 1.0E-9)
      {
        switch(TMath::Abs(candidate.PID))
        {
        case 11:
          fElectronOutputArray->emplace_back(candidate);
          break;
        case 13:
          fMuonOutputArray->emplace_back(candidate);
          break;
        default:
          fChargedHadronOutputArray->emplace_back(candidate);
        }
      }
      else
      {
        fNeutralOutputArray->emplace_back(candidate);
      }
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
      r_c = TMath::Hypot(x_c, y_c);

      // time of closest approach
      td = (phi_0 + TMath::ATan2(x_c, y_c)) / omega;

      // remove all the modulo pi that might have come from the atan
      pio = TMath::Abs(TMath::Pi() / omega);
      while(TMath::Abs(td) > 0.5 * pio)
      {
        td -= TMath::Sign(1.0, td) * pio;
      }

      vz = pz * c_light / e;

      // calculate coordinates of closest approach to z axis
      phid = phi_0 - omega * td;
      xd = x_c - r * TMath::Sin(phid);
      yd = y_c + r * TMath::Cos(phid);
      zd = z + vz * td;

      // momentum at closest approach
      px = pt * TMath::Cos(phid);
      py = pt * TMath::Sin(phid);

      const auto caParticleMomentum = ROOT::Math::PtEtaPhiEVector(pt, particleMomentum.Eta(), phid, particleMomentum.E());

      // calculate additional track parameters (correct for beamspot position)
      d0 = ((xd - bsx) * py - (yd - bsy) * px) / pt;
      dz = zd - bsz;
      ctgTheta = 1.0 / TMath::Tan(caParticleMomentum.Theta());

      // 3. time evaluation t = TMath::Min(t_r, t_z)
      //    t_r : time to exit from the sides
      //    t_z : time to exit from the front or the back
      t_z = (vz == 0.0) ? 1.0E99 : (TMath::Sign(fHalfLength, pz) - z) / vz;

      if(r_c + TMath::Abs(r) < fRadius)
      {
        // helix does not cross the cylinder sides
        t = t_z;
      }
      else
      {
        alpha = TMath::ACos((r * r + r_c * r_c - fRadius * fRadius) / (2 * TMath::Abs(r) * r_c));
        t_r = td + TMath::Abs(alpha / omega);

        t = TMath::Min(t_r, t_z);
      }

      // 4. position in terms of x(t), y(t), z(t)
      phi_t = phi_0 - omega * t;
      x_t = x_c - r * TMath::Sin(phi_t);
      y_t = y_c + r * TMath::Cos(phi_t);
      z_t = z + vz * t;
      r_t = TMath::Hypot(x_t, y_t);

      // lenght of the path from production to tracker
      l = t * TMath::Hypot(vz, r * omega);

      if(r_t > 0.0)
      {
        // store these variables before cloning
        if(particle == &candidate)
        {
          particle->D0 = d0 * 1.0E3;
          particle->DZ = dz * 1.0E3;
          particle->P = caParticleMomentum.P();
          particle->PT = pt;
          particle->CtgTheta = ctgTheta;
          particle->Phi = caParticleMomentum.Phi();
        }

        auto new_candidate = candidate;

        new_candidate.InitialPosition = particlePosition;
        new_candidate.Position.SetXYZT(x_t * 1.0E3, y_t * 1.0E3, z_t * 1.0E3, particlePosition.T() + t * c_light * 1.0E3);

        new_candidate.Momentum = caParticleMomentum;

        new_candidate.L = l * 1.0E3;

        new_candidate.Xd = xd * 1.0E3;
        new_candidate.Yd = yd * 1.0E3;
        new_candidate.Zd = zd * 1.0E3;

        new_candidate.AddCandidate(&candidate); // preserve parentage

        fOutputArray->emplace_back(new_candidate);
        switch(TMath::Abs(new_candidate.PID))
        {
        case 11:
          fElectronOutputArray->emplace_back(new_candidate);
          break;
        case 13:
          fMuonOutputArray->emplace_back(new_candidate);
          break;
        default:
          fChargedHadronOutputArray->emplace_back(new_candidate);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
