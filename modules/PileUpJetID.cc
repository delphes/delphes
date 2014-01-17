/** \class PileUpJetID
 *
 *  CMS PileUp Jet ID Variables, based on http://cds.cern.ch/record/1581583
 *
 *  \author S. Zenz, December 2013
 *
 *
 */

#include "modules/PileUpJetID.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

PileUpJetID::PileUpJetID() :
  fItJetInputArray(0),fTrackInputArray(0),fNeutralInputArray(0),fItVertexInputArray(0) 
{

}

//------------------------------------------------------------------------------

PileUpJetID::~PileUpJetID()
{

}

//------------------------------------------------------------------------------

void PileUpJetID::Init()
{
  fJetPTMin = GetDouble("JetPTMin", 20.0);
  fParameterR = GetDouble("ParameterR", 0.5);
  fUseConstituents = GetInt("UseConstituents", 0);

  fAverageEachTower = false; // for timing

  // import input array(s)

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "Calorimeter/eflowTracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  fNeutralInputArray = ImportArray(GetString("NeutralInputArray", "Calorimeter/eflowTowers"));
  fItNeutralInputArray = fNeutralInputArray->MakeIterator();
  
  fVertexInputArray = ImportArray(GetString("VertexInputArray", "PileUpMerger/vertices"));
  fItVertexInputArray = fVertexInputArray->MakeIterator();
  
  fZVertexResolution  = GetDouble("ZVertexResolution", 0.005)*1.0E3;
// create output array(s)

  fOutputArray = ExportArray(GetString("OutputArray", "jets"));

}

//------------------------------------------------------------------------------

void PileUpJetID::Finish()
{

  if(fItJetInputArray) delete fItJetInputArray;
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItNeutralInputArray) delete fItNeutralInputArray;
  if(fItVertexInputArray) delete fItVertexInputArray;

}

//------------------------------------------------------------------------------

void PileUpJetID::Process()
{
  Candidate *candidate, *constituent;
  TLorentzVector momentum, area;
  Double_t zvtx=0;

  Candidate *trk;

 // find z position of primary vertex
   
  fItVertexInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItVertexInputArray->Next())))
  {
    if(!candidate->IsPU)
    {
    zvtx = candidate->Position.Z();
    break;
    }
  }

  // loop over all input candidates
  fItJetInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItJetInputArray->Next())))
  {
    momentum = candidate->Momentum;
    area = candidate->Area;

    float sumpt = 0.;
    float sumptch = 0.;
    float sumptchpv = 0.;
    float sumptchpu = 0.;
    float sumdrsqptsq = 0.;
    float sumptsq = 0.;
    int nc = 0;
    int nn = 0;
    float pt_ann[5];

    for (int i = 0 ; i < 5 ; i++) {
      pt_ann[i] = 0.;
    }

    if (fUseConstituents) {
      TIter itConstituents(candidate->GetCandidates());
      while((constituent = static_cast<Candidate*>(itConstituents.Next()))) {
        float pt = constituent->Momentum.Pt();
        float dr = candidate->Momentum.DeltaR(constituent->Momentum);
        sumpt += pt;
        sumdrsqptsq += dr*dr*pt*pt;
        sumptsq += pt*pt;
        if (constituent->Charge == 0) {
	  // neutrals
	  nn++;
	} else {
	  // charged
	  if (constituent->IsPU && TMath::Abs(constituent->Position.Z()-zvtx) > fZVertexResolution) {
	    sumptchpu += pt;
	  } else {
	    sumptchpv += pt;
	  }
	  sumptch += pt;
	  nc++;
	}
	for (int i = 0 ; i < 5 ; i++) {
	  if (dr > 0.1*i && dr < 0.1*(i+1)) {
	    pt_ann[i] += pt;
	  }
	}
      }
    } else {
      // Not using constituents, using dr
      fItTrackInputArray->Reset();
       while ((trk = static_cast<Candidate*>(fItTrackInputArray->Next()))) {
	if (trk->Momentum.DeltaR(candidate->Momentum) < fParameterR) {
	  float pt = trk->Momentum.Pt();
	  sumpt += pt;
	  sumptch += pt;
	  if (trk->IsPU && TMath::Abs(trk->Position.Z()-zvtx) > fZVertexResolution) {
	    sumptchpu += pt;
	  } else {
	    sumptchpv += pt;
	  }
	  float dr = candidate->Momentum.DeltaR(trk->Momentum);
	  sumdrsqptsq += dr*dr*pt*pt;
	  sumptsq += pt*pt;
	  nc++;
	  for (int i = 0 ; i < 5 ; i++) {
	    if (dr > 0.1*i && dr < 0.1*(i+1)) {
              pt_ann[i] += pt;
	    }
	  }
	}
      }
      fItNeutralInputArray->Reset();
      while ((constituent = static_cast<Candidate*>(fItNeutralInputArray->Next()))) {
	if (constituent->Momentum.DeltaR(candidate->Momentum) < fParameterR) {
	  float pt = constituent->Momentum.Pt();
	  sumpt += pt;
	  float dr = candidate->Momentum.DeltaR(constituent->Momentum);
	  sumdrsqptsq += dr*dr*pt*pt;
	  sumptsq += pt*pt;
	  nn++;
	  for (int i = 0 ; i < 5 ; i++) {
	    if (dr > 0.1*i && dr < 0.1*(i+1)) {
	      pt_ann[i] += pt;
	    }
	  }
	}
      }
    }
          
    if (sumptch > 0.) {
      candidate->Beta = sumptchpu/sumptch;
      candidate->BetaStar = sumptchpv/sumptch;
    } else {
      candidate->Beta = -999.;
      candidate->BetaStar = -999.;
    }
    if (sumptsq > 0.) {
      candidate->MeanSqDeltaR = sumdrsqptsq/sumptsq;
    } else {
      candidate->MeanSqDeltaR = -999.;
    }
    candidate->NCharged = nc;
    candidate->NNeutrals = nn;
    if (sumpt > 0.) {
      candidate->PTD = TMath::Sqrt(sumptsq) / sumpt;
      for (int i = 0 ; i < 5 ; i++) {
        candidate->FracPt[i] = pt_ann[i]/sumpt;
      }
    } else {
      candidate->PTD = -999.;
      for (int i = 0 ; i < 5 ; i++) {
        candidate->FracPt[i] = -999.;
      }
    }

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------

