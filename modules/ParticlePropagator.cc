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

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>

using namespace std;

class ParticlePropagator: public DelphesModule
{
public:
  explicit ParticlePropagator(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fRadius(Steer<double>("Radius", 1.0)),
    fRadius2(fRadius * fRadius),
    fHalfLength(Steer<double>("HalfLength", 3.0)),
    fBz(Steer<double>("Bz", 0.0)),
    fRadiusMax(Steer<double>("RadiusMax", fRadius)),
    fHalfLengthMax(Steer<double>("HalfLengthMax", fHalfLength))
  {
    if(fRadius < 1.0E-2)
      throw std::runtime_error("ERROR: magnetic field radius is too low");
    if(fHalfLength < 1.0E-2)
      throw std::runtime_error("ERROR: magnetic field length is too low");
  }

  void Init() override
  {
    // import array with output from filter/classifier module
    fInputArray = ImportArray(Steer<std::string>("InputArray", "Delphes/stableParticles"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "stableParticles"));
    fNeutralOutputArray = ExportArray(Steer<std::string>("NeutralOutputArray", "neutralParticles"));
    fChargedHadronOutputArray = ExportArray(Steer<std::string>("ChargedHadronOutputArray", "chargedHadrons"));
    fElectronOutputArray = ExportArray(Steer<std::string>("ElectronOutputArray", "electrons"));
    fMuonOutputArray = ExportArray(Steer<std::string>("MuonOutputArray", "muons"));
    // import beamspot
    if(const std::string beamSpotArrayLabel = Steer<std::string>("BeamSpotInputArray", "BeamSpotFilter/beamSpotParticle");
      !beamSpotArrayLabel.empty() && GetFactory()->Has(beamSpotArrayLabel))
      fBeamSpotInputArray = ImportArray(beamSpotArrayLabel);
  }
  void Process() override;

private:
  const double fRadius;
  const double fRadius2;
  const double fHalfLength;
  const double fBz;
  const double fRadiusMax;
  const double fHalfLengthMax;

  CandidatesCollection fInputArray; //!

  CandidatesCollection fBeamSpotInputArray; //!

  CandidatesCollection fOutputArray; //!
  CandidatesCollection fNeutralOutputArray; //!
  CandidatesCollection fChargedHadronOutputArray; //!
  CandidatesCollection fElectronOutputArray; //!
  CandidatesCollection fMuonOutputArray; //!
};

//------------------------------------------------------------------------------

void ParticlePropagator::Process()
{
  fOutputArray->clear();
  fNeutralOutputArray->clear();
  fChargedHadronOutputArray->clear();
  fElectronOutputArray->clear();
  fMuonOutputArray->clear();

  const double c_light = 2.99792458E8;

  TLorentzVector beamSpotPosition;
  if(fBeamSpotInputArray && !fBeamSpotInputArray->empty())
  {
    Candidate &beamSpotCandidate = *((Candidate *)fBeamSpotInputArray->at(0));
    beamSpotPosition = beamSpotCandidate.Position;
  }

  Candidate *particle = nullptr;
  for(const Candidate *candidate : *fInputArray)
  {
    if(candidate->GetCandidates().empty())
      particle = const_cast<Candidate *>(candidate);
    else
      particle = static_cast<Candidate *>(candidate->GetCandidates().at(0));

    const TLorentzVector &particlePosition = particle->Position, &particleMomentum = particle->Momentum;

    const double x = particlePosition.X() * 1.0E-3, y = particlePosition.Y() * 1.0E-3, z = particlePosition.Z() * 1.0E-3;
    const double bsx = beamSpotPosition.X() * 1.0E-3, bsy = beamSpotPosition.Y() * 1.0E-3, bsz = beamSpotPosition.Z() * 1.0E-3;
    const double q = particle->Charge;

    // check that particle position is inside the cylinder
    if(std::hypot(x, y) > fRadiusMax || std::fabs(z) > fHalfLengthMax)
    {
      continue;
    }

    const double px = particleMomentum.Px(), py = particleMomentum.Py(), pz = particleMomentum.Pz();
    const double pt = particleMomentum.Pt();
    const double pt2 = particleMomentum.Perp2();
    const double e = particleMomentum.E();

    if(pt2 < 1.0E-9)
    {
      continue;
    }

    if(std::hypot(x, y) > fRadius || std::fabs(z) > fHalfLength)
    {
      Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());

      new_candidate->InitialPosition = particlePosition;
      new_candidate->Position = particlePosition;
      new_candidate->L = 0.0;

      new_candidate->Momentum = particleMomentum;
      new_candidate->AddCandidate(candidate);

      fOutputArray->emplace_back(new_candidate);
    }
    else if(std::fabs(q) < 1.0E-9 || std::fabs(fBz) < 1.0E-9)
    {
      // solve pt2*t^2 + 2*(px*x + py*y)*t - (fRadius2 - x*x - y*y) = 0
      const double tmp = px * y - py * x;
      const double t_r = (std::sqrt(pt2 * fRadius2 - tmp * tmp) - px * x - py * y) / pt2;
      const double t_z = (fHalfLength * pz / std::fabs(pz) - z) / pz;

      const double t = std::min(t_r, t_z);

      const double x_t = x + px * t, y_t = y + py * t, z_t = z + pz * t;

      const double l = std::sqrt((x_t - x) * (x_t - x) + (y_t - y) * (y_t - y) + (z_t - z) * (z_t - z));

      Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());

      new_candidate->InitialPosition = particlePosition;
      new_candidate->Position.SetXYZT(x_t * 1.0E3, y_t * 1.0E3, z_t * 1.0E3, particlePosition.T() + t * e * 1.0E3);
      new_candidate->L = l * 1.0E3;

      new_candidate->Momentum = particleMomentum;
      new_candidate->AddCandidate(candidate);

      fOutputArray->emplace_back(new_candidate);

      if(std::fabs(q) > 1.0E-9)
      {
        switch(std::abs(new_candidate->PID))
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
      else
      {
        fNeutralOutputArray->emplace_back(new_candidate);
      }
    }
    else
    {

      // 1. initial transverse momentum p_{T0}: Part->pt
      //    initial transverse momentum direction phi_0 = -atan(p_{X0} / p_{Y0})
      //    relativistic gamma: gamma = E / mc^2; gammam = gamma * m
      //    gyration frequency omega = q * Bz / (gammam)
      //    helix radius r = p_{T0} / (omega * gammam)

      const double gammam = e * 1.0E9 / (c_light * c_light); // gammam in [eV/c^2]
      const double omega = q * fBz / gammam; // omega is here in [89875518/s]
      const double r = pt / (q * fBz) * 1.0E9 / c_light; // in [m]

      const double phi_0 = std::atan2(py, px); // [rad] in [-pi, pi]

      // 2. helix axis coordinates
      const double x_c = x + r * std::sin(phi_0),
                   y_c = y - r * std::cos(phi_0),
                   r_c = std::hypot(x_c, y_c);

      // time of closest approach
      double td = (phi_0 + std::atan2(x_c, y_c)) / omega;

      // remove all the modulo pi that might have come from the atan
      double pio = std::fabs(M_PI / omega);
      while(std::fabs(td) > 0.5 * pio)
        td -= td / std::fabs(td) * pio;

      const double vz = pz * c_light / e;

      // calculate coordinates of closest approach to z axis
      const double phid = phi_0 - omega * td;
      const double xd = x_c - r * std::sin(phid),
                   yd = y_c + r * std::cos(phid),
                   zd = z + vz * td;

      // momentum at closest approach
      const double px = pt * std::cos(phid);
      const double py = pt * std::sin(phid);

      TLorentzVector particleMomentumAst;
      particleMomentumAst.SetPtEtaPhiE(pt, particleMomentum.Eta(), phid, particleMomentum.E());

      // calculate additional track parameters (correct for beamspot position)
      const double d0 = ((xd - bsx) * py - (yd - bsy) * px) / pt,
                   dz = zd - bsz;
      const double ctgTheta = 1.0 / std::tan(particleMomentumAst.Theta());

      // 3. time evaluation t = std::min(t_r, t_z)
      //    t_r : time to exit from the sides
      //    t_z : time to exit from the front or the back
      const double t_z = (vz == 0.0) ? 1.0E99 : (fHalfLength * pz / std::fabs(pz) - z) / vz;

      double t = 0.;
      if(r_c + std::fabs(r) < fRadius)
      {
        // helix does not cross the cylinder sides
        t = t_z;
      }
      else
      {
        const double alpha = std::acos((r * r + r_c * r_c - fRadius * fRadius) / (2 * std::fabs(r) * r_c));
        const double t_r = td + std::fabs(alpha / omega);

        t = std::min(t_r, t_z);
      }

      // 4. position in terms of x(t), y(t), z(t)
      const double phi_t = phi_0 - omega * t;
      const double x_t = x_c - r * std::sin(phi_t),
                   y_t = y_c + r * std::cos(phi_t),
                   z_t = z + vz * t,
                   r_t = std::hypot(x_t, y_t);

      // lenght of the path from production to tracker
      const double l = t * std::hypot(vz, r * omega);

      if(r_t > 0.0)
      {
        // store these variables before cloning
        if(particle == candidate)
        {
          particle->D0 = d0 * 1.0E3;
          particle->DZ = dz * 1.0E3;
          particle->P = particleMomentumAst.P();
          particle->PT = pt;
          particle->CtgTheta = ctgTheta;
          particle->Phi = particleMomentumAst.Phi();
        }

        Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());

        new_candidate->InitialPosition = particlePosition;
        new_candidate->Position.SetXYZT(x_t * 1.0E3, y_t * 1.0E3, z_t * 1.0E3, particlePosition.T() + t * c_light * 1.0E3);

        new_candidate->Momentum = particleMomentumAst;

        new_candidate->L = l * 1.0E3;

        new_candidate->Xd = xd * 1.0E3;
        new_candidate->Yd = yd * 1.0E3;
        new_candidate->Zd = zd * 1.0E3;

        new_candidate->AddCandidate(candidate);

        fOutputArray->emplace_back(new_candidate);
        switch(std::abs(new_candidate->PID))
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

REGISTER_MODULE("ParticlePropagator", ParticlePropagator);
