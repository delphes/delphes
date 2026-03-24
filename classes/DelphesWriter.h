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

#ifndef DelphesWriter_h
#define DelphesWriter_h

/** \class DelphesWriter
 *
 *  Base object to write event content file
 *
 *  \author L. Forthomme - AGH, Krakow
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"

class DelphesWriter: public DelphesModule
{
public:
  using DelphesModule::DelphesModule;
  virtual ~DelphesWriter() = default;

  void SetOutputFile(std::string_view outputFile) { fOutputFile = outputFile; }
  const std::string &GetOutputFile() const { return fOutputFile; }

  virtual void AddInfo(std::string_view name, double value) {}

  bool IsWriter() const override { return true; }

private:
  std::string fOutputFile;
};

#endif // DelphesWriter_h
