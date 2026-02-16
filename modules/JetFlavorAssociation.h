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

#ifndef JetFlavorAssociation_h
#define JetFlavorAssociation_h

/** \class JetFlavorAssociation
 *
 *  Find origin of jet and evaluate jet flavor
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class Candidate;
class DelphesFormula;

class ExRootSTLVectorFilter;
class PartonClassifier;
class ParticleLHEFClassifier;

class JetFlavorAssociation : public DelphesModule
{
public:
  JetFlavorAssociation();
  ~JetFlavorAssociation();

  void Init();
  void Process();
  void Finish();

  void GetAlgoFlavor(Candidate &jet, const std::vector<Candidate> &partonArray, const std::vector<Candidate> &partonLHEFArray);
  void GetPhysicsFlavor(Candidate &jet, const std::vector<Candidate> &partonArray, const std::vector<Candidate> &partonLHEFArray);

private:
  Double_t fDeltaR;

  PartonClassifier *fPartonClassifier; //!
  ParticleLHEFClassifier *fParticleLHEFClassifier; //!

  ExRootSTLVectorFilter *fPartonFilter;
  ExRootSTLVectorFilter *fParticleLHEFFilter;

  InputHandle<std::vector<Candidate> > fPartonInputArray; //!
  InputHandle<std::vector<Candidate> > fParticleInputArray; //!
  InputHandle<std::vector<Candidate> > fParticleLHEFInputArray; //!
  InputHandle<std::vector<Candidate> > fJetInputArray; //!

  ClassDef(JetFlavorAssociation, 1)
};

#endif
