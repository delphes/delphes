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

/** \class EICPIDDetector
 *
 *  Applies complex photon Id. Reconstructed photon candidtes are first separated into matched and non-matched to gen particles. 
 *  Non-matched pass the "fake" efficiency. Matched photons get further splitted into isolated and non-isolated (user can choose criterion for isolation)
 *  Isolated photons pass the "prompt" efficiency while the non-isolated pass the "non-prompt" efficiency
 *
 *  \author M. Selvaggi CERN
 *
 */

#include "modules/EICPIDDetector.h"

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
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"

#include "pid/barrelDIRC/src/barrelDirc.h"
#include "pid/quintRICH/src/CF4rich.h"
#include "pid/mRICH/src/mRICH.h"
#include "pid/tofBarrel/src/tofBarrel.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

EICPIDDetector::EICPIDDetector() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

EICPIDDetector::~EICPIDDetector()
{
}

//------------------------------------------------------------------------------

void EICPIDDetector::Init()
{

  // import input arrays
  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // PID Pair to be assessed
  ExRootConfParam param;
  Int_t size;
  param = GetParam("Hypotheses");
  size = param.GetSize();

  fHypo = static_cast<PID::type>(0);
  if (size == 2) {
    fPDG1 = abs(param[0].GetInt());
    fPDG2 = abs(param[1].GetInt());

    if (fPDG1 == 321) {
      if (fPDG2 == 211) {
    	fHypo = PID::pi_k;
      } else if (fPDG2 == 2212) {
    	fHypo = PID::k_p;
      } 
    }

  } else {
    // Bad parameter - do something intelligent here.
    std::cout << "Unable to retrieve Particle ID hypothesis pair." << std::endl;
  }


  fDetectorName = std::string(GetString("DetectorName", "barrelDirc"));

  // Common PID Detector parameters
  fTrackResolution = GetDouble("TrackResolution", 0.5); // mrad
  fTimeResolution = GetDouble("TimeResolution", 0.1); //ns
  fDetectorLength = GetDouble("DetectorLength", 1500); // mm
  fetaLow = GetDouble("EtaLow", -8.0);
  fetaHigh = GetDouble("EtaHigh", 8.0);

  // Barrel DIRC Parameters
  fQE = GetDouble("QuantumEfficiency", 0.0); // 0 = 27% for barrelDirc, 1 = 22% for barrelDirc

  // mRICH Parameters
  fPixelSize = GetDouble("PixelSize", 1.0); // 1.0 mm

  // CF4RICH Parameters
  

  // Build the detector object
  if (fDetectorName == "barrelDirc") {
    fPIDDetector = new barrelDirc(fTrackResolution,fTimeResolution,fQE,fetaLow,fetaHigh);
  } 
  else if (fDetectorName == "mRICH") {
    fPIDDetector = new mRICH(fTrackResolution,fTimeResolution, fPixelSize);
    //fPIDDetector = new mRICH(0.00175, 1, 3);
  }
  else if (fDetectorName == "CF4rich") {
    fPIDDetector = new CF4rich(fDetectorLength/10, fetaLow, fetaHigh, fPixelSize, fTrackResolution);
  }
  else if (fDetectorName == "tofBarrel") {
    fPIDDetector = new tofBarrel(100, fetaLow, fetaHigh, 10);
  } else {
    std::cout << "No valid EIC PID Detector technology was specified!" << std::endl;
    assert(1==0);
  }

  fPIDDetector->description();
  
  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void EICPIDDetector::Finish()
{
  if(fItInputArray) delete fItInputArray;
  if(fPIDDetector) delete fPIDDetector;
}

//------------------------------------------------------------------------------

void EICPIDDetector::Process()
{
  Candidate *candidate, *mother;
  Double_t pt, eta;
  Int_t true_id;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {

    mother = candidate;
    candidate = static_cast<Candidate *>(candidate->Clone());
    candidate->AddCandidate(mother);

    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = candidateMomentum.Eta();
    pt = candidateMomentum.Pt();
    true_id = candidate->PID;

    Float_t p = pt * TMath::CosH(eta);

    // Obtain the number of sigma separation for a given hypothesis pair for this track
    Bool_t valid = fPIDDetector->valid(eta, p);

    Double_t nsigma = -1.0;

    if (valid) {
      nsigma = fPIDDetector->numSigma(eta, p, fHypo);
      
      // Assume that Nsigma_Hypo1 = N_sigma_Hypo2, so that Nsigma_HypoX = Nsigma/Sqrt(2).
      nsigma = nsigma/TMath::Sqrt(2.0);
      //std::cout << std::scientific << "EICDetector nsigma = " << nsigma << std::fixed << std::endl;
    }
    


    int pid_reco = 0;
    int pid_true = TMath::Abs(candidate->PID);
    if (!valid || TMath::IsNaN(nsigma) || !TMath::Finite(nsigma)) {
      pid_reco = 0;
    } else {


      // Use accept/reject to assign a PID-detector identity to this track
      Double_t probability = 0.0;
      if (TMath::Abs(true_id) == fPDG1) {
	// We are selecting FOR this hypothesis, so use the core of a Gaussian as the probability
	probability = 1.0 - ROOT::Math::gaussian_pdf(nsigma);
      } else if (TMath::Abs(true_id) == fPDG2) {
	// We are trying to reject these using this detector, so the one-sided tail of the Gaussian probability applies
	probability = 1.0 - ROOT::Math::normal_cdf(nsigma);
      }
      
      // std::cout << "True ID = " << true_id << ", |eta| = " << TMath::Abs(eta) << ", p = " << pt*TMath::CosH(eta) 
      // 		<<", nsigma = " << std::scientific << nsigma << ", probability = " << std::fixed
      //  		<< std::scientific <<  probability << std::fixed << std::endl;
      
      // Create a PID value that is the concatenation of two 16-bit numbers.
      // The lowest 16 bits are the reconstructed PID
      // The highest 16 bits are the truth PID
      // Bitmasking and shifting can be used to get these separately.
      // For example, do the following to get the ... :
      //
      // * True PID: (Track.PID & 0xffff0000) >> 16) (Mask-select the highest 16 bits and shift right by 16 bits.
      // * Reco PID: (Track.PID & 0xffff) (Mask-select the lowest 16 bits)
      if (gRandom->Uniform(0, 1) < probability) {
	candidate = static_cast<Candidate *>(candidate->Clone());
	pid_reco = TMath::Abs(fPDG1);
	pid_true = TMath::Abs(candidate->PID);
      } else {
	pid_reco = TMath::Abs(fPDG2);
	pid_true = TMath::Abs(candidate->PID);
      }
    }
    int pid_all = pid_reco + (pid_true << 16);
    candidate->PID = pid_all;
    fOutputArray->Add(candidate);
  }
}
