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

#ifndef BoostedTagging_h
#define BoostedTagging_h

/** \class BoostedTagging
 *
 *  Determines the boosted origin of a large-R jet (W/Z/H/top),
 *  applies tagging efficiency (mis-identification rate) formulas
 *  and sets boosted-tagging flags
 *
 *  \author Sihyun Jeon - Boston University
 *
 */

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "classes/DelphesModule.h"

#include <map>

class TObjArray;
class DelphesFormula;

class ExRootFilter;
class BoostedTaggingClassifier;

class BoostedTagging: public DelphesModule
{
public:
  BoostedTagging();
  ~BoostedTagging();

  void Init();
  void Process();
  void Finish();

private:
  Int_t fBitNumber;

  Double_t fJetRadius;

  Double_t fDeltaR;

  Double_t fJetPTMin;

  Double_t fSoftDropMassMin;

  Double_t fSoftDropMassMax;

#if !defined(__CINT__) && !defined(__CLING__)
  std::map<Int_t, DelphesFormula *> fEfficiencyMap; //!
#endif

  BoostedTaggingClassifier *fClassifier = nullptr; //!

  ExRootFilter *fFilter = nullptr;

  TIterator *fItJetInputArray = nullptr; //!

  const TObjArray *fParticleInputArray = nullptr; //!

  const TObjArray *fJetInputArray = nullptr; //!

  ClassDef(BoostedTagging, 1)
};

//------------------------------------------------------------------------------

class BoostedTaggingClassifier: public ExRootClassifier
{
public:
  BoostedTaggingClassifier(const TObjArray *array);

  Int_t GetCategory(TObject *object);

  Double_t fEtaMax, fPTMin;

  const TObjArray *fParticleInputArray = nullptr;
};

#endif
