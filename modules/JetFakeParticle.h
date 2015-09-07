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

#ifndef JetFakeParticle_h
#define JetFakeParticle_h


/** \class JetFakeParticle
 *
 *  Converts jet into particle with some PID,
 *  according to parametrized probability.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;
class DelphesFormula;

class JetFakeParticle: public DelphesModule
{
public:

  JetFakeParticle();
  ~JetFakeParticle();

  void Init();
  void Process();
  void Finish();

private:

  #if !defined(__CINT__) && !defined(__CLING__)
  typedef std::map< Int_t, DelphesFormula * > TFakeMap; //!
  TFakeMap fEfficiencyMap;
  #endif

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fElectronOutputArray; //!
  TObjArray *fMuonOutputArray; //!
  TObjArray *fPhotonOutputArray; //!
  TObjArray *fJetOutputArray; //!

  ClassDef(JetFakeParticle, 1)
};

#endif
