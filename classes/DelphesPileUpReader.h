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

#ifndef DelphesPileUpReader_h
#define DelphesPileUpReader_h

/** \class DelphesPileUpReader
 *
 *  Reads pile-up binary file
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include <array>
#include <memory>
#include <string_view>

#include <stdint.h>
#include <stdio.h>

class DelphesXDRReader;

class DelphesPileUpReader
{
public:
  explicit DelphesPileUpReader(std::string_view fileName);
  ~DelphesPileUpReader();

  bool ReadParticle(int32_t &pid,
    float &x, float &y, float &z, float &t,
    float &px, float &py, float &pz, float &e);

  bool ReadEntry(int64_t entry);

  int64_t GetEntries() const { return fEntries; }

private:
  static constexpr size_t kIndexSize = 10000000;
  static constexpr size_t kBufferSize = 1000000;
  static constexpr size_t kRecordSize = 9;

  int64_t fEntries{0};

  int32_t fEntrySize{0};
  int32_t fCounter{0};

  FILE *fPileUpFile{nullptr};
  std::array<uint8_t, kIndexSize * 8> fIndex;
  std::array<uint8_t, kBufferSize * kRecordSize * 4> fBuffer;

  const std::unique_ptr<DelphesXDRReader> fInputReader;
  const std::unique_ptr<DelphesXDRReader> fIndexReader;
  const std::unique_ptr<DelphesXDRReader> fBufferReader;
};

#endif // DelphesPileUpReader_h
