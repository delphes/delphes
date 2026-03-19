/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2026  Universite catholique de Louvain (UCL), Belgium
 *                           AGH University of Krakow, Poland
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

#ifndef DelphesReader_h
#define DelphesReader_h

/** \class DelphesReader
 *
 *  Base object to read event content file
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve, L. Forthomme - AGH, Krakow
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModuleFactory.h"

#include <TStopwatch.h>

class DelphesFactory;

class ExRootProgressBar;
class ExRootTreeBranch;

class DelphesReader
{
public:
  DelphesReader() = default;
  virtual ~DelphesReader();

  virtual void LoadInputFile(std::string_view inputFile) {}
  virtual void Clear() = 0;
  virtual bool EventReady() = 0;

  void SetSkipEvents(long long);
  void SetMaxEvents(long long);
  unsigned long long EventCounter() const { return fEventCounter; }

  virtual bool ReadEvent(DelphesFactory *factory,
    CandidatesCollection &allParticleOutputArray,
    CandidatesCollection &stableParticleOutputArray,
    CandidatesCollection &partonOutputArray);

  virtual void AnalyzeEvent(ExRootTreeBranch *branch, TStopwatch *procStopWatch) = 0;
  virtual void AnalyzeWeight(ExRootTreeBranch *branch) {}

protected:
  virtual bool ReadBlock(DelphesFactory *factory,
    CandidatesCollection &allParticleOutputArray,
    CandidatesCollection &stableParticleOutputArray,
    CandidatesCollection &partonOutputArray) { return false; }

  FILE *fInputFile{nullptr};
  std::unique_ptr<ExRootProgressBar> fProgressBar;
  TStopwatch fReadStopWatch;

  long long fEventCounter{0ll};
  long long fMaxEvents{-1ll};
  long long fSkipEvents{-1ll};
};

/// Add an event reader to the list of handled modules
#define REGISTER_READER(name, obj)                                                 \
  struct BUILDER_NAME(obj)                                                         \
  {                                                                                \
    BUILDER_NAME(obj)() { DelphesReaderFactory::Get().RegisterModule<obj>(name); } \
  };                                                                               \
  static const BUILDER_NAME(obj) gDelphesReader##obj;                              \
  static_assert(true, "")

/// A documentation generator factory
DEFINE_FACTORY(DelphesReaderFactory, DelphesReader, "Event readers factory");

#endif // DelphesReader_h
