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

#ifndef DelphesLHEFReader_h
#define DelphesLHEFReader_h

/** \class DelphesLHEFReader
 *
 *  Reads LHEF file
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"

#include <stdio.h>
#include <utility>
#include <vector>

#include "classes/DelphesReader.h"

class TDatabasePDG;

class DelphesLHEFReader: public DelphesReader
{
public:
  DelphesLHEFReader();
  ~DelphesLHEFReader();

  void SetInputFile(FILE *inputFile);

  void Clear() override;
  bool EventReady() override;

  bool ReadBlock(DelphesFactory *factory,
    CandidatesCollection &allParticleOutputArray,
    CandidatesCollection &stableParticleOutputArray,
    CandidatesCollection &partonOutputArray) override;

  void AnalyzeEvent(ExRootTreeBranch *branch, long long eventNumber,
    TStopwatch *readStopWatch, TStopwatch *procStopWatch) override;

  void AnalyzeWeight(ExRootTreeBranch *branch);

private:
  void AnalyzeParticle(DelphesFactory *factory,
    CandidatesCollection &allParticleOutputArray,
    CandidatesCollection &stableParticleOutputArray,
    CandidatesCollection &partonOutputArray);

  FILE *fInputFile;

  char *fBuffer;

  TDatabasePDG *fPDG;

  bool fEventReady;

  int fEventCounter;

  int fParticleCounter, fProcessID;
  double fCrossSection, fWeight, fScalePDF, fAlphaQCD, fAlphaQED;

  int fPID, fStatus, fM1, fM2, fC1, fC2;
  double fPx, fPy, fPz, fE, fMass;

  std::vector<std::pair<int, double> > fWeightList;
};

#endif // DelphesLHEFReader_h
