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

#ifndef ParticlePropagator_h
#define ParticlePropagator_h

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

#include "classes/DelphesModel.h"
#include "classes/DelphesModule.h"

class Candidate;

class ParticlePropagator : public DelphesModule
{
public:
  ParticlePropagator() = default;

  void Init();
  void Process();
  void Finish();

private:
  Double_t fRadius, fRadius2, fRadiusMax, fHalfLength, fHalfLengthMax;
  Double_t fBz;

  InputHandle<std::vector<Candidate> > fInputArray; //!
  InputHandle<std::vector<Candidate> > fBeamSpotInputArray; //!
  OutputHandle<std::vector<Candidate> > fOutputArray; //!
  OutputHandle<std::vector<Candidate> > fNeutralOutputArray; //!
  OutputHandle<std::vector<Candidate> > fChargedHadronOutputArray; //!
  OutputHandle<std::vector<Candidate> > fElectronOutputArray; //!
  OutputHandle<std::vector<Candidate> > fMuonOutputArray; //!

  ClassDef(ParticlePropagator, 1)
};

#endif
