/** \class BeamSpotFilter
 *
 *  Extracts beam spot
 *
 *  \author Michele Selvaggi
 *
 */

#include "modules/BeamSpotFilter.h"

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

BeamSpotFilter::BeamSpotFilter() :
  fItInputArray(0), fVertexArray(0), fItVertexArray(0)
{
}

//------------------------------------------------------------------------------

BeamSpotFilter::~BeamSpotFilter()
{
}

//------------------------------------------------------------------------------

void BeamSpotFilter::Init()
{

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "filteredParticles"));

  // determine mode -- do we define the beamspot using truth-level info (Candidate::IsPU),
  // or do we define it from the vertex with the largest sum(pt^2) of charged particles?
  fMode = GetInt("PrimaryVertexMode",0);

  // resolution in meters for vertexing -- only relevant for fMode==1 use case
  fResolution = GetDouble("Resolution", 1E-06);
}

//------------------------------------------------------------------------------

void BeamSpotFilter::Finish()
{
  if(fItInputArray) delete fItInputArray;
  if(fVertexArray) delete fVertexArray;
  if(fItVertexArray) delete fItVertexArray;
}

//------------------------------------------------------------------------------

void BeamSpotFilter::Process()
{
  Candidate *candidate;
  Bool_t passed = false;

  fItInputArray->Reset();

  if(fMode == 1)
  {
    // produce the vertex objects
    ProduceVertices();
    fItVertexArray->Reset();
    Double_t sumpt2 = 0.;
    Double_t sumpt2_max = 0.;
    Int_t index = -1;
    Int_t counter = 0;

    while((candidate = static_cast<Candidate *>(fItVertexArray->Next())))
    {
      sumpt2 = candidate->GenSumPT2;
      if(sumpt2 > sumpt2_max)
      {
        sumpt2_max = sumpt2;
        index = counter;
      }
      counter++;
    }

    candidate = static_cast<Candidate *>(fInputArray->At(index));
    fOutputArray->Add(candidate);
  }
  else // the "standard" functionality, using pileup info
    {

    while((candidate = static_cast<Candidate *>(fItInputArray->Next())) && !passed)
    {
      if(candidate->IsPU == 0) passed = true;
      fOutputArray->Add(candidate);
    }
}
}

void BeamSpotFilter::ProduceVertices()
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
     fItVertexArray = fVertexArray->MakeIterator();
     fItVertexArray->Reset();
     while((vertex = static_cast<Candidate *>(fItVertexArray->Next())))
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
        fVertexArray->Add(vertex);
        nvtx++;
      }
   }
}