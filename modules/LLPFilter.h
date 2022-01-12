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

//------------------------------------------------------------------------------

#ifndef LLPFilter_h
#define LLPFilter_h


/** \class LLPFilter
 *
 *  Filter LLPs with particular PDG ID/status and calculate the EM and hadronic energy of LLP based on decay particles
 *  The classification of EM and hadronic energy of LLP is based on instructions from the HEPData entry for the CMS paper searching
 *  for neutral LLPs in the CMS endcap muon detectors: https://www.hepdata.net/record/104408
 *
 *  \author Christina Wang
 *
 */

#include "classes/DelphesModule.h"
#include <vector>

class TIterator;
class TObjArray;

class LLPFilter: public DelphesModule
{
public:
  LLPFilter();
  ~LLPFilter();

  void Init();
  void Process();
  void Finish();

private:
  Double_t fPTMin; //!
  Bool_t fRequireDecayRegion;
  Double_t fDecayRegionRMax;
  Double_t fDecayRegionRMin;
  Double_t fDecayRegionZMax;
  Double_t fDecayRegionZMin;
  Double_t fDecayRegionEtaMax;
  Double_t fDecayRegionEtaMin;
  Int_t fDaughterNumber;
  Bool_t fInvert; //!
  Bool_t fRequireStatus; //!
  Int_t fStatus; //!
  Bool_t fRequireCharge; //!
  Int_t fCharge; //!
  Bool_t fRequireNotPileup; //!

  std::vector<Int_t> fPdgCodes;

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TIterator *fItParticleInputArray;
  const TObjArray *fParticleInputArray;

  TObjArray *fOutputArray; //!

  ClassDef(LLPFilter, 1)
};

#endif
