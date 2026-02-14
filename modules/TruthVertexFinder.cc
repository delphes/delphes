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

/** \class TruthVertexFinder
 *
 *  Merges particles from pile-up sample into event
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/TruthVertexFinder.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesPileUpReader.h"
#include "classes/DelphesTF2.h"

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

void TruthVertexFinder::Init()
{
  fResolution = GetDouble("Resolution", 1E-06); // resolution in meters
  // import input array
  ImportArray(GetString("InputArray", "Delphes/stableParticles"), fInputArray);
  // create output array
  ExportArray(fVertexOutputArray, GetString("VertexOutputArray", "vertices"));
}

//------------------------------------------------------------------------------

void TruthVertexFinder::Finish()
{
}

//------------------------------------------------------------------------------

void TruthVertexFinder::Process()
{
  Int_t nvtx = -1;
  Float_t pt;

  ROOT::Math::XYZTVector vertexPosition(0., 0., 0., 0.);

  nvtx = 0;
  for(const auto &candidate : *fInputArray)
  {

    const auto &candidatePosition = candidate.Position;
    const auto &candidateMomentum = candidate.Momentum;

    pt = candidateMomentum.Pt();

    // check whether vertex already included, if so add particle
    Bool_t old_vertex = false;
    for(auto &vertex : *fVertexOutputArray)
    {
      const auto &vertexPosition = vertex.Position;
      // check whether spatial difference is < 1 um, in that case assume it is the same vertex
      if(TMath::Abs((candidatePosition.P() - vertexPosition.P())) < fResolution * 1.E3)
      {
        old_vertex = true;
        vertex.AddCandidate(const_cast<Candidate *>(&candidate)); //TODO: ensure const-qualification
        if(TMath::Abs(candidate.Charge) > 0)
        {
          vertex.ClusterNDF += 1;
          vertex.GenSumPT2 += pt * pt;
        }
      }
    }

    // else fill new vertex
    if(!old_vertex)
    {
      auto *vertex = GetFactory()->NewCandidate();
      vertex->Position = candidatePosition;
      vertex->ClusterIndex = nvtx;

      if(TMath::Abs(candidate.Charge) > 0)
      {
        vertex->ClusterNDF = 1;
        vertex->GenSumPT2 = pt * pt;
      }
      else
      {
        vertex->ClusterNDF = 0;
        vertex->GenSumPT2 = 0.;
      }
      fVertexOutputArray->emplace_back(*vertex);
      nvtx++;
    }
  }
}

//------------------------------------------------------------------------------
