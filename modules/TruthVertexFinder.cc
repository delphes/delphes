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

class TruthVertexFinder: public DelphesModule
{
public:
  explicit TruthVertexFinder(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fResolution(Steer<double>("Resolution", 1E-06)) // resolution in meters
  {
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "Delphes/stableParticles"));
    fVertexOutputArray = ExportArray(Steer<std::string>("VertexOutputArray", "vertices"));
  }
  void Process() override;

private:
  const double fResolution; //!

  CandidatesCollection fInputArray; //!
  CandidatesCollection fVertexOutputArray; //!
};

//------------------------------------------------------------------------------

void TruthVertexFinder::Process()
{
  fVertexOutputArray->clear();

  DelphesFactory *factory = GetFactory();

  int nvtx = 0;
  for(Candidate *const &candidate : *fInputArray)
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;

    const double pt = candidateMomentum.Pt();

    // check whether vertex already included, if so add particle
    bool old_vertex = false;
    for(Candidate *const &vertex : *fVertexOutputArray)
    {
      const TLorentzVector &vertexPosition = vertex->Position;
      // check whether spatial difference is < 1 um, in that case assume it is the same vertex
      if(std::fabs((candidatePosition.P() - vertexPosition.P())) < fResolution * 1.E3)
      {
        old_vertex = true;
        vertex->AddCandidate(candidate);
        if(std::fabs(candidate->Charge) > 0)
        {
          vertex->ClusterNDF += 1;
          vertex->GenSumPT2 += pt * pt;
        }
      }
    }

    // else fill new vertex
    if(!old_vertex)
    {
      Candidate *vertex = factory->NewCandidate();
      vertex->Position = candidatePosition;
      vertex->ClusterIndex = nvtx;

      if(std::fabs(candidate->Charge) > 0)
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
