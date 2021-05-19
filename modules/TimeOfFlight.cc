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

TimeOfFlight::TimeOfFlight() :
  fItTrackInputArray(0), fItVertexInputArray(0)
{
}

//------------------------------------------------------------------------------

TimeOfFlight::~TimeOfFlight()
{
}

//------------------------------------------------------------------------------

void TimeOfFlight::Init()
{

  // method to compute vertex time
  fVertexTimeMode = GetInt("VertexTimeMode", 0);

  // import track input array
  fTrackInputArray = ImportArray(GetString("TrackInputArray", "MuonMomentumSmearing/muons"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  // import vertex input array
  fVertexInputArray = ImportArray(GetString("VertexInputArray", "TruthVertexFinder/vertices"));
  fItVertexInputArray = fVertexInputArray->MakeIterator();

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void TimeOfFlight::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItVertexInputArray) delete fItVertexInputArray;
}

//------------------------------------------------------------------------------

void TimeOfFlight::Process()
{
  Candidate *candidate, *particle, *vertex, *constituent, *mother;
  Double_t ti, t_truth, tf;
  Double_t l, tof, beta, p,  mass;

  const Double_t c_light = 2.99792458E8;

  // first compute momenta of vertices based on reconstructed tracks
  ComputeVertexMomenta();

  fItTrackInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItTrackInputArray->Next())))
  {

    particle = static_cast<Candidate *>(candidate->GetCandidates()->At(0));

    const TLorentzVector &candidateInitialPosition = particle->Position;
    const TLorentzVector &candidateInitialPositionSmeared = candidate->InitialPosition;
    const TLorentzVector &candidateFinalPosition = candidate->Position;
    const TLorentzVector &candidateMomentum = particle->Momentum;

    // time at vertex from MC truth
    t_truth = candidateInitialPosition.T() * 1.0E-3 / c_light;

    if (candidate->Position.Vect().Mag() < 5.) continue;

    // various options on how to calculate the vertex time
    ti=0;
    switch (fVertexTimeMode)
	  {
    	case 0:
    	{
        // assume ti from MC truth
        // most aggressive, we are cheating and assume we can perfectly reconstruct time of primary and secondary vertices
        ti = t_truth;
        break;
    	}
    	case 1:
    	{
        // always assume t=0, most conservative assumption
        // reasonable assumption for particles originating from PV, if beamSpot has small time spread compared to timing resolution
        // probably bad assumption for particles from highly displaced vertices (i.e Ks)
        ti=0;
        break;
    	}
    	case 2:
    	{
        // same as 2 but attempt at estimate beta from vertex mass and momentum
        beta = 1.;
        fItVertexInputArray->Reset();
        while((vertex = static_cast<Candidate *>(fItVertexInputArray->Next())))
        {
          TIter itGenParts(vertex->GetCandidates());
          itGenParts.Reset();

          while((constituent = static_cast<Candidate *>(itGenParts.Next())))
          {
            if (particle == constituent)
            {
              beta = vertex->Momentum.Beta();
              break;
            }
          }
        } // end vertex  loop

        // track displacement to be possibily replaced by vertex fitted position
        ti = candidateInitialPositionSmeared.Vect().Mag() * 1.0E-3 /(beta*c_light);
      }
      break;
  	}

    p = candidateMomentum.P();

    // this quantity has already been smeared by another module
    tf = candidateFinalPosition.T() * 1.0E-3 / c_light;

    // calculate time-of-flight
    tof = tf - ti;
    // path length of the full helix
    l = candidate->L * 1.0E-3;

    // particle velocity
    beta = l/(c_light*tof);

    // calculate particle mass (i.e particle ID)
    mass = 0.;
    if (beta<1) mass = p* TMath::Sqrt(1/(beta*beta) - 1);

    mother    = candidate;
    candidate = static_cast<Candidate*>(candidate->Clone());

    // update time at vertex based on option
    candidate->InitialPosition.SetT(ti * 1.0E3 * c_light);

    // update particle mass based on TOF-based PID (commented for now, assume this calculation is done offline)
    //candidate->Momentum.SetVectM(candidateMomentum.Vect(), mass);

    candidate->AddCandidate(mother);
    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------

void TimeOfFlight::ComputeVertexMomenta()
{
  Candidate *track, *constituent, *particle, *vertex;

  fItVertexInputArray->Reset();
  while((vertex = static_cast<Candidate *>(fItVertexInputArray->Next())))
  {

    TIter itGenParts(vertex->GetCandidates());
    itGenParts.Reset();

    while((constituent = static_cast<Candidate *>(itGenParts.Next())))
    {
      fItTrackInputArray->Reset();
      while((track = static_cast<Candidate *>(fItTrackInputArray->Next())))
      {
        // get gen part that generated track
        particle = static_cast<Candidate *>(track->GetCandidates()->At(0));
        if (particle == constituent)
        {
          vertex->Momentum += track->Momentum;
        }

      } // end track loop
    } // end vertex consitutent loop
  } // end vertex  loop
}
