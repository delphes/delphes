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
#include "classes/DelphesFactory.h"
#include "classes/DelphesPileUpReader.h"
#include "classes/DelphesTF2.h"

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

TruthVertexFinder::TruthVertexFinder() :
fItInputArray(0), fItOutputArray(0)
{
}

//------------------------------------------------------------------------------

TruthVertexFinder::~TruthVertexFinder()
{
}

//------------------------------------------------------------------------------

void TruthVertexFinder::Init()
{

  fResolution = GetDouble("Resolution", 1E-06); // resolution in meters
  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output arrays
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));
  //fItOutputArray = fVertexOutputArray->MakeIterator();
}

//------------------------------------------------------------------------------

void TruthVertexFinder::Finish()
{
  if(fItInputArray) delete fItInputArray;
  if(fItOutputArray) delete fItOutputArray;
}

//------------------------------------------------------------------------------

void TruthVertexFinder::Process()
{
  Int_t nvtx = -1;
  Float_t pt;
  Candidate *candidate, *vertex;
  DelphesFactory *factory;


  fItInputArray->Reset();

  factory = GetFactory();
  vertex = factory->NewCandidate();

  TLorentzVector vertexPosition(0., 0., 0., 0.);

  nvtx=0;
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {

     const TLorentzVector &candidatePosition = candidate->Position;
     const TLorentzVector &candidateMomentum = candidate->Momentum;

     pt = candidateMomentum.Pt();
     
     // check whether vertex already included, if so add particle
     Bool_t old_vertex=false;
     fItOutputArray = fVertexOutputArray->MakeIterator();
     fItOutputArray->Reset();
     while((vertex = static_cast<Candidate *>(fItOutputArray->Next())))
     {
        const TLorentzVector &vertexPosition = vertex->Position;
        // check whether spatial difference is < 1 um, in that case assume it is the same vertex
        if ( TMath::Abs((candidatePosition.P() - vertexPosition.P())) < fResolution*1.E3)
        {
           old_vertex=true;
           vertex->AddCandidate(candidate);
           if (TMath::Abs(candidate->Charge) > 0)
           {
              vertex->ClusterNDF += 1;
              vertex->GenSumPT2 += pt*pt;
           }
        }
     }

     // else fill new vertex
     if (!old_vertex)
     {
        vertex = factory->NewCandidate();
        vertex->Position = candidatePosition;
        vertex->ClusterIndex = nvtx;

        if (TMath::Abs(candidate->Charge) > 0)
        {
           vertex->ClusterNDF = 1;
           vertex->GenSumPT2 = pt*pt;
        }
        else
        {
           vertex->ClusterNDF = 0;
           vertex->GenSumPT2 = 0.;
        }
        fVertexOutputArray->Add(vertex);
        nvtx++;
      }
   }
}

//------------------------------------------------------------------------------
