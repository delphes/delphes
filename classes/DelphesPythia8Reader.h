/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2026  Universite catholique de Louvain (UCLouvain), Belgium
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

#ifndef DelphesPythia8Reader_h
#define DelphesPythia8Reader_h

/** \class DelphesPythia8Reader
 *
 *  Reads Pythia8 data
 *
 *  \author P. Demin - UCLouvain, Louvain-la-Neuve
 *
 */

#include <string>

#include <cstdio>

#include "TObject.h"

class TObjArray;
class TStopwatch;
class TDatabasePDG;
class ExRootTreeBranch;
class DelphesFactory;

namespace Pythia8
{
class Pythia;
class CombineMatchingInput;
} // namespace Pythia8

class DelphesPythia8Reader: public TObject
{
public:
  DelphesPythia8Reader();
  ~DelphesPythia8Reader();

  void OpenInputFile(const char *inputFileName);
  void CloseInputFile();

  void SetInputFile(FILE *inputFile);

  bool ReadLHEF();
  std::string FileNameLHEF();

  void Clear(Option_t *option = "") override;
  bool EventReady();

  int EventNumber();

  bool ReadEvent(DelphesFactory *factory,
    TObjArray *allParticleOutputArray,
    TObjArray *stableParticleOutputArray,
    TObjArray *partonOutputArray);

  void AnalyzeEvent(ExRootTreeBranch *branch, long long eventNumber,
    TStopwatch *readStopWatch = 0, TStopwatch *procStopWatch = 0);

  void AnalyzeWeight(ExRootTreeBranch *branch);

private:
  void AnalyzeParticles(DelphesFactory *factory,
    TObjArray *allParticleOutputArray,
    TObjArray *stableParticleOutputArray,
    TObjArray *partonOutputArray);

  Pythia8::Pythia *fPythia = nullptr;
  Pythia8::CombineMatchingInput *fCombine = nullptr;

  TDatabasePDG *fPDG = nullptr;

  Long64_t fEventCounter = 0;
  Long64_t fErrorCounter = 0;

  Long64_t fNumberOfEvents = 0;
  Long64_t fTimesAllowErrors = 0;

  Int_t fFrameType = 0;
  std::string fFileNameLHEF;

  Bool_t fSpareFlag1 = false;
  Int_t fSpareMode1 = 0;
  Double_t fSpareParm1 = 0.0;
  Double_t fSpareParm2 = 0.0;

  bool fEventReady = false;

  ClassDef(DelphesPythia8Reader, 1)
};

#endif // DelphesPythia8Reader_h
