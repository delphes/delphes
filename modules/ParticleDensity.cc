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

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"

#include <TH2F.h>

using namespace std;

class ParticleDensity: public DelphesModule
{
public:
  ParticleDensity() = default;

  void Init() override;
  void Process() override;

private:
  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!

  Bool_t fUseMomentumVector; // !
  std::unique_ptr<TH2F> fHisto; //!
};

//------------------------------------------------------------------------------

void ParticleDensity::Init()
{
  // import input array
  fInputArray = ImportArray(GetString("InputArray", "FastJetFinder/jets"));

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));

  // create multiplicity histogram
  ExRootConfParam paramEta = GetParam("EtaBins");
  const Long_t sizeEta = paramEta.GetSize();
  Int_t nbinsEta = sizeEta - 1;
  std::vector<Float_t> binsEta(static_cast<size_t>(sizeEta));
  for(Int_t i = 0; i < sizeEta; ++i)
  {
    binsEta[i] = paramEta[i].GetDouble();
  }

  ExRootConfParam paramPhi = GetParam("PhiBins");
  const Long_t sizePhi = paramPhi.GetSize();
  Int_t nbinsPhi = sizePhi - 1;
  std::vector<Float_t> binsPhi(static_cast<size_t>(sizePhi));
  for(Int_t i = 0; i < sizePhi; ++i)
  {
    binsPhi[i] = paramPhi[i].GetDouble();
  }

  fHisto = std::make_unique<TH2F>("hParticleDensity", ";#eta;#varphi;d^{2}N/d#etad#varphi", nbinsEta, binsEta.data(), nbinsPhi, binsPhi.data());

  fUseMomentumVector = GetBool("UseMomentumVector", false);
}

//------------------------------------------------------------------------------

void ParticleDensity::Process()
{
  fOutputArray->clear();

  fHisto->Reset();

  // loop over all input candidates to fill histogram
  for(Candidate *const &candidate : *fInputArray)
  {
    if(fUseMomentumVector)
      fHisto->Fill(candidate->Momentum.Eta(), candidate->Momentum.Phi());
    else
      fHisto->Fill(candidate->Position.Eta(), candidate->Position.Phi());
  }

  // normalise by bin width
  fHisto->Scale(1., "width");

  // loop over all input candidates to assign multiplicity
  for(Candidate *const &candidate : *fInputArray)
  {
    Int_t ieta, iphi;
    if(fUseMomentumVector)
    {
      ieta = fHisto->GetXaxis()->FindBin(candidate->Momentum.Eta());
      iphi = fHisto->GetYaxis()->FindBin(candidate->Momentum.Phi());
    }
    else
    {
      ieta = fHisto->GetXaxis()->FindBin(candidate->Position.Eta());
      iphi = fHisto->GetYaxis()->FindBin(candidate->Position.Phi());
    }
    candidate->ParticleDensity = fHisto->GetBinContent(ieta, iphi);
    fOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("ParticleDensity", ParticleDensity);
