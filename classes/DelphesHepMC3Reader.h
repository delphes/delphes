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
  DelphesHepMC3Reader();

  void LoadInputFile(std::string_view) override;
  void Clear() override;
  bool EventReady() override;

  void AnalyzeEvent(ExRootTreeBranch *branch, TStopwatch *procStopWatch) override;
  void AnalyzeWeight(ExRootTreeBranch *branch) override;

private:
  bool ReadBlock(DelphesFactory *factory,
    CandidatesCollection &allParticleOutputArray,
    CandidatesCollection &stableParticleOutputArray,
    CandidatesCollection &partonOutputArray) override;

  void AnalyzeVertex(DelphesFactory *factory, int code, Candidate *candidate = 0);

  void AnalyzeParticle(DelphesFactory *factory);

  void FinalizeParticles(CandidatesCollection &allParticleOutputArray,
    CandidatesCollection &stableParticleOutputArray,
    CandidatesCollection &partonOutputArray);

  FILE *fInputFile{nullptr};

  std::array<char, 16384> fBuffer;

  TDatabasePDG *fPDG;

  int fEventNumber, fMPI, fProcessID, fVertexCounter, fParticleCounter;
  double fScale, fAlphaQCD, fAlphaQED;

  double fMomentumCoefficient, fPositionCoefficient;

  std::vector<double> fWeights;

  double fCrossSection, fCrossSectionError;

  int fID1, fID2;
  double fX1, fX2, fScalePDF, fPDF1, fPDF2;

  int fVertexCode, fVertexStatus;
  double fX, fY, fZ, fT;

  int fParticleCode, fPID, fParticleStatus, fOutVertexCode;
  double fPx, fPy, fPz, fE, fMass;

  std::vector<std::pair<std::shared_ptr<TLorentzVector>, CandidatesCollection> > fVertices;
  std::vector<int> fParticles;

  std::map<int, int> fInVertexMap;
  std::map<int, int> fOutVertexMap;

  std::map<int, std::pair<int, int> > fMotherMap;
  std::map<int, std::pair<int, int> > fDaughterMap;
};

#endif // DelphesHepMC3Reader_h
