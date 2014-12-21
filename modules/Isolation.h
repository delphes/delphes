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

#ifndef Isolation_h
#define Isolation_h

/** \class Isolation
 *
 *  Sums transverse momenta of isolation objects (tracks, calorimeter towers, etc)
 *  within a DeltaR cone around a candidate and calculates fraction of this sum
 *  to the candidate's transverse momentum. outputs candidates that have
 *  the transverse momenta fraction within (PTRatioMin, PTRatioMax].
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TObjArray;

class ExRootFilter;
class IsolationClassifier;

class Isolation: public DelphesModule
{
public:

  Isolation();
  ~Isolation();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fDeltaRMax;

  Double_t fPTRatioMax;

  Double_t fPTSumMax;

  Bool_t fUsePTSum;

  IsolationClassifier *fClassifier; //!

  ExRootFilter *fFilter;

  TIterator *fItIsolationInputArray; //!

  TIterator *fItCandidateInputArray; //!

  TIterator *fItRhoInputArray; //!

  const TObjArray *fIsolationInputArray; //!

  const TObjArray *fCandidateInputArray; //!

  const TObjArray *fRhoInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(Isolation, 1)
};

#endif
