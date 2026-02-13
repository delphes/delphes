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

/** \class ParticleDensity
 *
 *  This module calculates the particle multiplicity density in eta-phi bins.
 *  It then assigns the value to the candidates according to the candidate eta.
 *
 *  \author R. Preghenella - INFN, Bologna
 *
 */

#include "modules/ParticleDensity.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TFormula.h"
#include "TH2F.h"
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

void ParticleDensity::Init()
{
  // import input array
  GetFactory()->EventModel()->Attach(GetString("InputArray", "FastJetFinder/jets"), fInputArray); // I/O
  // create output array
  GetFactory()->EventModel()->Book(fOutputArray, GetString("OutputArray", "tracks"));

  // create multiplicity histogram

  ExRootConfParam paramEta = GetParam("EtaBins");
  const Long_t sizeEta = paramEta.GetSize();
  Int_t nbinsEta = sizeEta - 1;
  Float_t binsEta[sizeEta];
  for(Int_t i = 0; i < sizeEta; ++i)
  {
    binsEta[i] = paramEta[i].GetDouble();
  }

  ExRootConfParam paramPhi = GetParam("PhiBins");
  const Long_t sizePhi = paramPhi.GetSize();
  Int_t nbinsPhi = sizePhi - 1;
  Float_t binsPhi[sizePhi];
  for(Int_t i = 0; i < sizePhi; ++i)
  {
    binsPhi[i] = paramPhi[i].GetDouble();
  }

  fHisto = new TH2F("hParticleDensity", ";#eta;#varphi;d^{2}N/d#etad#varphi", nbinsEta, binsEta, nbinsPhi, binsPhi);

  fUseMomentumVector = GetBool("UseMomentumVector", false);
}

//------------------------------------------------------------------------------

void ParticleDensity::Finish()
{
  if(fHisto) delete fHisto;
}

//------------------------------------------------------------------------------

void ParticleDensity::Process()
{
  fHisto->Reset();

  // loop over all input candidates to fill histogram
  for(const auto &candidate : *fInputArray)
  {
    if(fUseMomentumVector)
      fHisto->Fill(candidate.Momentum.Eta(), candidate.Momentum.Phi());
    else
      fHisto->Fill(candidate.Position.Eta(), candidate.Position.Phi());
  }

  // normalise by bin width
  fHisto->Scale(1., "width");

  // loop over all input candidates to assign multiplicity
  for(auto &candidate : *fInputArray)
  {
    Int_t ieta, iphi;
    if(fUseMomentumVector)
    {
      ieta = fHisto->GetXaxis()->FindBin(candidate.Momentum.Eta());
      iphi = fHisto->GetYaxis()->FindBin(candidate.Momentum.Phi());
    }
    else
    {
      ieta = fHisto->GetXaxis()->FindBin(candidate.Position.Eta());
      iphi = fHisto->GetYaxis()->FindBin(candidate.Position.Phi());
    }
    candidate.ParticleDensity = fHisto->GetBinContent(ieta, iphi);
    fOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------
