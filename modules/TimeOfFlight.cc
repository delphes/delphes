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

/** \class TimeOfFlight
  *
  *  Calculates Time-Of-Flight
  *
  *  \author Michele Selvaggi - CERN
  *
 */

#include "modules/TimeOfFlight.h"

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

void TimeOfFlight::Init()
{

  // method to compute vertex time
  fVertexTimeMode = GetInt("VertexTimeMode", 0);

  // import track input array
  GetFactory()->EventModel()->Attach(GetString("InputArray", "MuonMomentumSmearing/muons"), fInputArray);

  // import vertex input array
  GetFactory()->EventModel()->Attach(GetString("VertexInputArray", "TruthVertexFinder/vertices"), fVertexInputArray);

  // create output array
  ExportArray(fOutputArray, GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void TimeOfFlight::Finish()
{
}

//------------------------------------------------------------------------------

void TimeOfFlight::Process()
{
  Candidate *constituent;
  Double_t ti, t_truth, tf;
  Double_t l, tof, beta;

  const Double_t c_light = 2.99792458E8;

  // first compute momenta of vertices based on reconstructed tracks
  ComputeVertexMomenta();

  for(const auto &candidate : *fInputArray) //TODO: ensure const-qualification of consumers
  {

    auto *particle = static_cast<Candidate *>(const_cast<Candidate &>(candidate).GetCandidates()->At(0));

    const auto &candidateInitialPosition = particle->Position;
    const auto &candidateInitialPositionSmeared = candidate.InitialPosition;
    const auto &candidateFinalPosition = candidate.Position;

    // time at vertex from MC truth
    t_truth = candidateInitialPosition.T() * 1.0E-3 / c_light;

    // various options on how to calculate the vertex time
    ti = 0;
    switch(fVertexTimeMode)
    {
    case 0: {
      // assume ti from MC truth
      // most aggressive, we are cheating and assume we can perfectly reconstruct time of primary and secondary vertices
      ti = t_truth;
      break;
    }
    case 1: {
      // always assume t=0, most conservative assumption
      // reasonable assumption for particles originating from PV, if beamSpot has small time spread compared to timing resolution
      // probably bad assumption for particles from highly displaced vertices (i.e Ks)
      ti = 0;
      break;
    }
    case 2: {
      // same as 2 but attempt at estimate beta from vertex mass and momentum
      beta = 1.;
      for(auto &vertex : *fVertexInputArray) //TODO: ensure const-qualification of consumers
      {
        TIter itGenParts(vertex.GetCandidates());
        itGenParts.Reset();

        while((constituent = static_cast<Candidate *>(itGenParts.Next())))
        {
          if(particle == constituent)
          {
            beta = vertex.Momentum.Beta();
            break;
          }
        }
      } // end vertex  loop

      // track displacement to be possibily replaced by vertex fitted position
      ti = std::sqrt(candidateInitialPositionSmeared.Vect().Mag2()) * 1.0E-3 / (beta * c_light);
    }
    break;
    }

    // this quantity has already been smeared by another module
    tf = candidateFinalPosition.T() * 1.0E-3 / c_light;

    // calculate time-of-flight
    tof = tf - ti;
    // path length of the full helix
    l = candidate.L * 1.0E-3;

    // particle velocity
    beta = l / (c_light * tof);

    // calculate particle mass (i.e particle ID)
    auto new_candidate = candidate;

    // update time at vertex based on option
    new_candidate.InitialPosition.SetE(ti * 1.0E3 * c_light);

    // update particle mass based on TOF-based PID (commented for now, assume this calculation is done offline)
    //new_candidate.Momentum.SetVectM(candidateMomentum.Vect(), mass);

    new_candidate.AddCandidate(const_cast<Candidate *>(&candidate)); // keep parentage
    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------

void TimeOfFlight::ComputeVertexMomenta()
{
  Candidate *constituent;

  for(auto &vertex : *fVertexInputArray)
  {
    TIter itGenParts(vertex.GetCandidates());
    itGenParts.Reset();

    while((constituent = static_cast<Candidate *>(itGenParts.Next())))
    {
      for(auto &track : *fInputArray)
      {
        // get gen part that generated track
        auto *particle = static_cast<Candidate *>(track.GetCandidates()->At(0));
        if(particle == constituent)
          vertex.Momentum += track.Momentum;
      } // end track loop
    } // end vertex consitutent loop
  } // end vertex  loop
}
