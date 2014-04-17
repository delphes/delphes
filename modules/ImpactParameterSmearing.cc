/** \class ImpactParameterSmearing
 *
 *  Performs transverse impact parameter smearing.
 *
 *  $Date: 2014-16-03 14:57:44 +0100   
 *
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */
 

#include "modules/ImpactParameterSmearing.h"

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

ImpactParameterSmearing::ImpactParameterSmearing() :
  fFormula(0), fItInputArray(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

ImpactParameterSmearing::~ImpactParameterSmearing()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void ImpactParameterSmearing::Init()
{
  // read resolution formula

  fFormula->Compile(GetString("ResolutionFormula", "0.0"));

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "TrackMerger/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void ImpactParameterSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void ImpactParameterSmearing::Process()
{
  Candidate *candidate, *particle, *mother;
  Double_t xd, yd, zd, dxy, dz, sx, sy, sz, ddxy, ddz;
  Double_t pt, eta, px, py, ang_mom;

 // cout<<"New event"<<endl;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
  
    //take momentum before smearing (otherwise apply double smearing on dxy)
    particle = static_cast<Candidate*>(candidate->GetCandidates()->At(0));
  
    const TLorentzVector &candidateMomentum = particle->Momentum;
    //  const TLorentzVector &candidateMomentum = candidate->Momentum;
    
    eta = candidateMomentum.Eta();
    pt = candidateMomentum.Pt();
    px = candidateMomentum.Px();
    py = candidateMomentum.Py();
      
    // calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
    xd =  candidate->Xd;
    yd =  candidate->Yd;
    zd =  candidate->Zd;
   
    // calculate smeared values   
    sx = gRandom->Gaus(0,fFormula->Eval(pt, eta));
    sy = gRandom->Gaus(0,fFormula->Eval(pt, eta));
    sz = gRandom->Gaus(0,fFormula->Eval(pt, eta));
     
    xd += sx;
    yd += sy;
    zd += sz; 
     
    // calculate impact paramater (after-smearing)
    ang_mom = (xd*py - yd*px);
    dxy = ang_mom/pt;
    dz = zd;
  
    ddxy = gRandom->Gaus(0,fFormula->Eval(pt, eta));
    ddz = gRandom->Gaus(0,fFormula->Eval(pt, eta));
  
    //fill smeared values in candidate
    mother = candidate;
    
    candidate = static_cast<Candidate*>(candidate->Clone());
    candidate->Xd = xd;
    candidate->Yd = yd;
    candidate->Zd = zd;
    
    candidate->Dxy = dxy;
    candidate->SDxy = ddxy;
      
    candidate->AddCandidate(mother);
    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
