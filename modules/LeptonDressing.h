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

#ifndef LeptonDressing_h
#define LeptonDressing_h

/** \class LeptonDressing
 *
 *  \author P. Demin && A. Mertens - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;

class LeptonDressing: public DelphesModule
{
public:

  LeptonDressing();
  ~LeptonDressing();
  
  void Init();
  void Process();
  void Finish();

private:

  Double_t fDeltaR;
  
  TIterator *fItDressingInputArray; //!
  
  TIterator *fItCandidateInputArray; //!

  const TObjArray *fDressingInputArray; //!
  
  const TObjArray *fCandidateInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(LeptonDressing, 1)
};

#endif
