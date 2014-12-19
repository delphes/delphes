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

#ifndef TauTagging_h
#define TauTagging_h

/** \class TauTagging
 *
 *  Determines origin of jet,
 *  applies b-tagging efficiency (miss identification rate) formulas
 *  and sets b-tagging flags 
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <map>

class TObjArray;
class DelphesFormula;

class ExRootFilter;
class TauTaggingPartonClassifier;

class TauTagging: public DelphesModule
{
public:

  TauTagging();
  ~TauTagging();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fDeltaR;

#if !defined(__CINT__) && !defined(__CLING__)
  std::map< Int_t, DelphesFormula * > fEfficiencyMap; //!
#endif
  
  TauTaggingPartonClassifier *fClassifier; //!
  
  ExRootFilter *fFilter;

  TIterator *fItPartonInputArray; //!
  
  TIterator *fItJetInputArray; //!

  const TObjArray *fParticleInputArray; //!

  const TObjArray *fPartonInputArray; //!
  
  const TObjArray *fJetInputArray; //!

  ClassDef(TauTagging, 1)
};

#endif
