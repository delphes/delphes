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

#ifndef TaggingParticlesSkimmer_h
#define TaggingParticlesSkimmer_h

/** \class TaggingParticlesSkimmer
 *
 *  Filters particle collection by only keeping gen particles useful for b and tau tagging.
    tau particles are replaced by their "visible decay".
 *
 *  \author M. Selvaggi
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;

class ExRootFilter;
class TauTaggingPartonClassifier;

class TaggingParticlesSkimmer: public DelphesModule
{
public:

  TaggingParticlesSkimmer();
  ~TaggingParticlesSkimmer();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fPTMin; //!
  Double_t fEtaMax; //!
  
  TauTaggingPartonClassifier *fClassifier; //!
  
  ExRootFilter *fFilter;

  TIterator *fItPartonInputArray; //!
 
  const TObjArray *fPartonInputArray; //!
  const TObjArray *fParticleInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(TaggingParticlesSkimmer, 1)
};

#endif
