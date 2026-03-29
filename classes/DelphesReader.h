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
#include "classes/DelphesModule.h"
#include "classes/DelphesModuleFactory.h"

#include <TStopwatch.h>

class DelphesFactory;

class ExRootProgressBar;

class DelphesReader: public DelphesModule
{
public:
  explicit DelphesReader(const DelphesParameters &readerParams = DelphesParameters{}) : DelphesModule(readerParams) {}
  virtual ~DelphesReader();

  virtual void LoadInputFile(std::string_view inputFile) {}
  virtual void Reset() {};
  virtual void Clear() = 0;

  void SetFactory(DelphesFactory *) override;

  void SetSkipEvents(long long);
  void SetMaxEvents(long long);
  unsigned long long EventCounter() const { return fEventCounter; }

  virtual bool ReadEvent();
  virtual void AnalyzeEvent(TStopwatch * /*procStopWatch*/) {}

  bool IsReader() const override { return true; }

protected:
  virtual bool EventReady() { return false; }
  virtual bool ReadBlock() { return false; }

  CandidatesCollection fAllParticleOutputArray;
  CandidatesCollection fStableParticleOutputArray;
  CandidatesCollection fPartonOutputArray;

  FILE *fInputFile{nullptr};
  std::unique_ptr<ExRootProgressBar> fProgressBar;
  TStopwatch fReadStopWatch;

  unsigned long long fEventCounter{0ll};
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
DEFINE_FACTORY(DelphesReaderFactory, DelphesReader, "event readers");

#endif // DelphesReader_h
