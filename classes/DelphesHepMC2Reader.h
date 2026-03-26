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

#ifndef DelphesHepMC2Reader_h
#define DelphesHepMC2Reader_h

/** \class DelphesHepMC2Reader
 *
 *  Reads HepMC file
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include <map>
#include <vector>

#include <stdio.h>

#include "classes/DelphesReader.h"

class TDatabasePDG;

class HepMCEvent;
class Weight;

class DelphesHepMC2Reader: public DelphesReader
{
public:
  explicit DelphesHepMC2Reader(const DelphesParameters &);

  void SetFactory(DelphesFactory *factory) override;

  void LoadInputFile(std::string_view) override;
  void Clear() override;
  bool EventReady() override;

  void AnalyzeEvent(TStopwatch *procStopWatch) override;

private:
  bool ReadBlock() override;

  void AnalyzeParticle();

  void FinalizeParticles();

  std::shared_ptr<HepMCEvent> fEventObject;
  std::shared_ptr<std::vector<Weight> > fWeights;

  FILE *fInputFile{nullptr};
  std::array<char, 16384> fBuffer;

  TDatabasePDG *fPDG{nullptr};

  int fVertexCounter{-1};
  double fMomentumCoefficient{1.}, fPositionCoefficient{1.};

  std::vector<int> fState;

  int fOutVertexCode{-1}, fVertexID{-1}, fInCounter{-1}, fOutCounter{-1};
  double fX{0.}, fY{0.}, fZ{0.}, fT{0.};

  int fParticleCode{0}, fPID{0}, fStatus{0}, fInVertexCode{-1};
  double fPx{0.}, fPy{0.}, fPz{0.}, fE{0.}, fMass{0.}, fTheta{0.}, fPhi{0.};

  int fParticleCounter{0};

  std::map<int, std::pair<int, int> > fMotherMap;
  std::map<int, std::pair<int, int> > fDaughterMap;
};

#endif // DelphesHepMC2Reader_h
