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

/** \class TrackCountingBTagging
 *
 *  b-tagging algorithm based on counting tracks with large impact parameter
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>

using namespace std;

class TrackCountingBTagging: public DelphesModule
{
public:
  explicit TrackCountingBTagging(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fBitNumber(Steer<int>("BitNumber", 0)),
    //
    fPtMin(Steer<double>("TrackPtMin", 1.0)),
    fDeltaR(Steer<double>("DeltaR", 0.3)),
    fIPmax(Steer<double>("TrackIPMax", 2.0)),
    //
    fSigMin(Steer<double>("SigMin", 6.5)),
    fNtracks(Steer<int>("Ntracks", 3)),
    //
    fUse3D(Steer<bool>("Use3D", false))
  {
  }

  void Init() override
  {
    fTrackInputArray = ImportArray(Steer<std::string>("TrackInputArray", "Calorimeter/eflowTracks"));
    fJetInputArray = ImportArray(Steer<std::string>("JetInputArray", "FastJetFinder/jets"));
  }
  void Process() override;

private:
  const int fBitNumber;

  const double fPtMin;
  const double fDeltaR;
  const double fIPmax;
  const double fSigMin;
  const int fNtracks;
  const bool fUse3D;

  CandidatesCollection fTrackInputArray; //!
  CandidatesCollection fJetInputArray; //!
};

//------------------------------------------------------------------------------

void TrackCountingBTagging::Process()
{
  // loop over all input jets
  for(Candidate *const &jet : *fJetInputArray)
  {
    const TLorentzVector &jetMomentum = jet->Momentum;
    const double jpx = jetMomentum.Px(), jpy = jetMomentum.Py(), jpz = jetMomentum.Pz();

    // loop over all input tracks
    int count = 0;
    for(Candidate *const &track : *fTrackInputArray)
    {
      if(count >= fNtracks) break; // stop once we have enough tracks
      const TLorentzVector &trkMomentum = track->Momentum;
      const double tpt = trkMomentum.Pt();
      if(tpt < fPtMin) continue;

      const double d0 = std::fabs(track->D0);
      if(d0 > fIPmax) continue;

      const double dr = jetMomentum.DeltaR(trkMomentum);
      if(dr > fDeltaR) continue;

      const double xd = track->Xd,
                   yd = track->Yd,
                   zd = track->Zd;
      const double dd0 = std::fabs(track->ErrorD0);
      const double dz = std::fabs(track->DZ);
      const double ddz = std::fabs(track->ErrorDZ);

      double sip = 0.;
      int sign = -1;
      if(fUse3D)
      {
        sign = (jpx * xd + jpy * yd + jpz * zd > 0.0) ? 1 : -1;
        //add transverse and longitudinal significances in quadrature
        sip = sign * std::sqrt(std::pow(d0 / dd0, 2) + std::pow(dz / ddz, 2));
      }
      else
      {
        sign = (jpx * xd + jpy * yd > 0.0) ? 1 : -1;
        sip = sign * d0 / std::fabs(dd0);
      }

      if(sip > fSigMin) count++;
    }

    // set BTag flag to true if count >= Ntracks
    jet->BTag |= (count >= fNtracks) << fBitNumber;
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("TrackCountingBTagging", TrackCountingBTagging);
