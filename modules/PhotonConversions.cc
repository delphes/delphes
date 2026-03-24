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
 *  Converts photons into e+ e- pairs according to material ditribution in the detector.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesCylindricalFormula.h"
#include "classes/DelphesModule.h"

#include <TF1.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TVector3.h>

using namespace std;

class PhotonConversions: public DelphesModule
{
public:
  explicit PhotonConversions(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fRadius(Steer<double>("Radius", 1.0)),
    fRadius2(fRadius * fRadius),
    fHalfLength(Steer<double>("HalfLength", 3.0)),
    fEtaMin(Steer<double>("EtaMin", 2.0)),
    fEtaMax(Steer<double>("EtaMax", 5.0)),
    fStep(Steer<double>("Step", 0.1)), // in meters
    fConversionMap(std::make_unique<DelphesCylindricalFormula>()),
    fDecayXsec(std::make_unique<TF1>("decayXsec", "1.0 - 4.0/3.0 * x * (1.0 - x)", 0.0, 1.0))
  {
    fConversionMap->Compile(Steer<std::string>("ConversionMap", "0.0"));
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "Delphes/stableParticles"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "stableParticles"));
  }
  void Process() override;

private:
  const double fRadius;
  const double fRadius2;
  const double fHalfLength;
  const double fEtaMin;
  const double fEtaMax;
  const double fStep;

  const std::unique_ptr<DelphesCylindricalFormula> fConversionMap; //!
  const std::unique_ptr<TF1> fDecayXsec; //!

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void PhotonConversions::Process()
{
  fOutputArray->clear();

  for(Candidate *const &candidate : *fInputArray)
  {
    if(candidate->PID != 22)
      fOutputArray->emplace_back(candidate);
    else
    {
      const TLorentzVector &candidatePosition = candidate->Position,
                           &candidateMomentum = candidate->Momentum;
      const double x = candidatePosition.X() * 1.0E-3,
                   y = candidatePosition.Y() * 1.0E-3,
                   z = candidatePosition.Z() * 1.0E-3;

      // check that particle position is inside the cylinder
      if(std::hypot(x, y) > fRadius || std::fabs(z) > fHalfLength) continue;

      const double px = candidateMomentum.Px(),
                   py = candidateMomentum.Py(),
                   pz = candidateMomentum.Pz(),
                   pt = candidateMomentum.Pt(),
                   pt2 = candidateMomentum.Perp2(),
                   eta = candidateMomentum.Eta(),
                   phi = candidateMomentum.Phi(),
                   e = candidateMomentum.E();

      if(eta < fEtaMin || eta > fEtaMax) continue;

      // solve pt2*t^2 + 2*(px*x + py*y)*t - (fRadius2 - x*x - y*y) = 0
      double tmp = px * y - py * x;
      const double discr2 = pt2 * fRadius2 - tmp * tmp;
      if(discr2 < 0.) // no solutions
        continue;

      tmp = px * x + py * y;
      const double discr = std::sqrt(discr2);
      const double t1 = (-tmp + discr) / pt2, t2 = (-tmp - discr) / pt2;
      double t = (t1 < 0.0) ? t2 : t1;

      if(const double z_t = z + pz * t; std::fabs(z_t) > fHalfLength)
      {
        const double t3 = (+fHalfLength - z) / pz, t4 = (-fHalfLength - z) / pz;
        t = (t3 < 0.0) ? t4 : t3;
      }

      // final position
      const double x_t = x + px * t, y_t = y + py * t, z_t = z + pz * t,
                   r_t = std::hypot(x_t, y_t, z_t);

      // here starts conversion code
      const int nsteps = int(r_t / fStep);

      const double dt = t / nsteps;

      bool converted = false;

      double x_i = x, y_i = y, z_i = z;
      for(int i = 0; i < nsteps; ++i)
      {
        x_i += px * dt;
        y_i += py * dt;
        z_i += pz * dt;
        TVector3 pos_i(x_i, y_i, z_i);

        // convert photon position into cylindrical coordinates, cylindrical r,phi,z !!

        const double r_i = std::hypot(x_i, y_i), phi_i = pos_i.Phi();

        // read conversion rate/meter from card
        const double rate = fConversionMap->Eval(r_i, phi_i, z_i);

        // convert into conversion probability
        const double p_conv = 1 - std::exp(-7.0 / 9.0 * fStep * rate);

        // case conversion occurs
        if(gRandom->Uniform() < p_conv)
        {
          converted = true;

          // generate x1 and x2, the fraction of the photon energy taken resp. by e+ and e-
          const double x1 = fDecayXsec->GetRandom(), x2 = 1 - x1;

          Candidate *ep = static_cast<Candidate *>(candidate->Clone());
          Candidate *em = static_cast<Candidate *>(candidate->Clone());

          ep->Position.SetXYZT(x_i * 1.0E3, y_i * 1.0E3, z_i * 1.0E3, candidatePosition.T() + nsteps * dt * e * 1.0E3);
          em->Position.SetXYZT(x_i * 1.0E3, y_i * 1.0E3, z_i * 1.0E3, candidatePosition.T() + nsteps * dt * e * 1.0E3);

          ep->Momentum.SetPtEtaPhiE(x1 * pt, eta, phi, x1 * e);
          em->Momentum.SetPtEtaPhiE(x2 * pt, eta, phi, x2 * e);

          ep->PID = -11;
          em->PID = 11;

          ep->Charge = 1.0;
          em->Charge = -1.0;

          ep->IsFromConversion = 1;
          em->IsFromConversion = 1;

          fOutputArray->emplace_back(em);
          fOutputArray->emplace_back(ep);

          break;
        }
      }
      if(!converted) fOutputArray->emplace_back(candidate);
    }
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("PhotonConversions", PhotonConversions);
