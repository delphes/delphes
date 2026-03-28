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

#ifndef DelphesSTDHEPReader_h
#define DelphesSTDHEPReader_h

/** \class DelphesSTDHEPReader
 *
 *  Reads STDHEP file
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include <stdint.h>
#include <stdio.h>

#include "classes/DelphesReader.h"
#include "classes/DelphesXDRReader.h"

class TDatabasePDG;

class DelphesSTDHEPReader: public DelphesReader
{
public:
  enum STDHEPBlock
  {
    GENERIC = 0,
    FILEHEADER = 1,
    EVENTTABLE = 2,
    EVENTHEADER = 4,
    MCFIO_STDHEP = 101,
    MCFIO_STDHEPBEG = 106,
    MCFIO_STDHEPEND = 107,
    MCFIO_STDHEP4 = 201
  };

  explicit DelphesSTDHEPReader(const DelphesParameters &);

  void SetFactory(DelphesFactory *factory) override;

  void LoadInputFile(std::string_view) override;
  void Reset() override;
  void Clear() override;
  bool EventReady() override;

  void AnalyzeEvent(TStopwatch *procStopWatch) override;

private:
  bool ReadBlock() override;

  void AnalyzeParticles();

  void SkipBytes(int size);
  void SkipArray(int elsize);

  void ReadFileHeader();
  void ReadEventTable();
  void ReadEventHeader();
  void ReadSTDCM1();
  void ReadSTDHEP();
  void ReadSTDHEP4();

  FILE *fInputFile{nullptr};

  std::shared_ptr<LHEFEvent> fEventInfo;

  DelphesXDRReader fReader[7];

  std::array<uint8_t, 1000000 * 96 + 24> fBuffer;

  TDatabasePDG *fPDG{nullptr};

  int32_t fBlockType{-1}, fEventNumber, fEventSize;
  std::vector<double> fScale;
};

#endif // DelphesSTDHEPReader_h
