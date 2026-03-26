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

class DelphesWriter;

class DelphesLHEFReader: public DelphesReader
{
public:
  explicit DelphesLHEFReader(const DelphesParameters &);

  void SetFactory(DelphesFactory *factory) override;

  void LoadInputFile(std::string_view) override;
  void Clear() override;
  bool EventReady() override { return fEventReady; }

  void AnalyzeEvent(TStopwatch *procStopWatch) override;

private:
  bool ReadBlock() override;
  void AnalyzeParticle();

  std::shared_ptr<LHEFEvent> fEventInfo;
  std::shared_ptr<std::vector<LHEFWeight> > fWeightInfo;

  std::array<char, 16384> fBuffer;

  TDatabasePDG *fPDG{nullptr};

  bool fEventReady{false};

  int fParticleCounter{-1};

  int fPID{0}, fStatus{0}, fM1{-1}, fM2{-1}, fC1{-1}, fC2{-1};
  double fPx{0.}, fPy{0.}, fPz{0.}, fE{0.}, fMass{0.};
};

#endif // DelphesLHEFReader_h
