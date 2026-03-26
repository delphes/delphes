/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2021  Universite catholique de Louvain (UCL), Belgium
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

#ifndef DelphesHepMC3Reader_h
#define DelphesHepMC3Reader_h

/** \class DelphesHepMC3Reader
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
class TLorentzVector;

class Candidate;

class DelphesHepMC3Reader: public DelphesReader
{
public:
  explicit DelphesHepMC3Reader(const DelphesParameters &);

  void SetFactory(DelphesFactory *factory) override;

  void LoadInputFile(std::string_view) override;
  void Clear() override;
  bool EventReady() override;

  void AnalyzeEvent(TStopwatch *procStopWatch) override;

private:
  bool ReadBlock() override;

  void AnalyzeVertex(int code, Candidate *candidate = 0);

  void AnalyzeParticle();

  void FinalizeParticles();

  std::shared_ptr<HepMCEvent> fEventObject;
  std::shared_ptr<std::vector<Weight> > fWeightsObject;

  FILE *fInputFile{nullptr};

  std::array<char, 16384> fBuffer;

  TDatabasePDG *fPDG{nullptr};

  double fMomentumCoefficient{1.}, fPositionCoefficient{1.};

  int fVertexCounter{0}, fParticleCounter{0};
  std::vector<double> fWeights;

  double fX{0.}, fY{0.}, fZ{0.}, fT{0.};

  int fParticleCode{0}, fPID{0}, fParticleStatus{0}, fOutVertexCode{0};
  double fPx{0.}, fPy{0.}, fPz{0.}, fE{0.}, fMass{0.};

  std::vector<std::pair<std::shared_ptr<TLorentzVector>, CandidatesCollection> > fVertices;
  std::vector<int> fParticles;

  std::map<int, int> fInVertexMap;
  std::map<int, int> fOutVertexMap;

  std::map<int, std::pair<int, int> > fMotherMap;
  std::map<int, std::pair<int, int> > fDaughterMap;
};

#endif // DelphesHepMC3Reader_h
