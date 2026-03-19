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
 *  Produces list of MC truth vertices
 *
 *  \author M. Selvaggi - CERN
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TMath.h>

using namespace std;

class TruthVertexFinder: public DelphesModule
{
public:
  TruthVertexFinder() = default;

  void Init() override;
  void Process() override;

private:
  Double_t fResolution; //!

  CandidatesCollection fInputArray; //!
  CandidatesCollection fVertexOutputArray; //!
};

//------------------------------------------------------------------------------

void TruthVertexFinder::Init()
{
  fResolution = GetDouble("Resolution", 1E-06); // resolution in meters

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/stableParticles"));

  // create output array
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));
}

//------------------------------------------------------------------------------

void TruthVertexFinder::Process()
{
  fVertexOutputArray->clear();

  Int_t nvtx = -1;
  Float_t pt;
  DelphesFactory *factory = nullptr;

  factory = GetFactory();

  TLorentzVector vertexPosition(0., 0., 0., 0.);

  nvtx = 0;
  for(const auto &candidate : *fInputArray)
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;

    pt = candidateMomentum.Pt();

    // check whether vertex already included, if so add particle
    Bool_t old_vertex = false;
    for(auto &vertex : *fVertexOutputArray)
    {
      const TLorentzVector &vertexPosition = vertex->Position;
      // check whether spatial difference is < 1 um, in that case assume it is the same vertex
      if(TMath::Abs((candidatePosition.P() - vertexPosition.P())) < fResolution * 1.E3)
      {
        old_vertex = true;
        vertex->AddCandidate(candidate);
        if(TMath::Abs(candidate->Charge) > 0)
        {
          vertex->ClusterNDF += 1;
          vertex->GenSumPT2 += pt * pt;
        }
      }
    }

    // else fill new vertex
    if(!old_vertex)
    {
      auto *vertex = factory->NewCandidate();
      vertex->Position = candidatePosition;
      vertex->ClusterIndex = nvtx;

      if(TMath::Abs(candidate->Charge) > 0)
      {
        vertex->ClusterNDF = 1;
        vertex->GenSumPT2 = pt * pt;
      }
      else
      {
        vertex->ClusterNDF = 0;
        vertex->GenSumPT2 = 0.;
      }
      fVertexOutputArray->emplace_back(vertex);
      nvtx++;
    }
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("TruthVertexFinder", TruthVertexFinder);
