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

/** \class PyDelphesConfig
 *
 *  Python configuration holder
 *
 *  \author L. Forthomme - AGH, Kraków
 *
 */

#ifndef DelphesPython_PyDelphesConfReader_h
#define DelphesPython_PyDelphesConfReader_h

//#include <pybind11/embed.h>
#include <pybind11/pybind11.h>

#include "classes/DelphesConfReader.h"
#include "classes/DelphesParameters.h"

class PyDelphesConfReader: public DelphesConfReader
{
public:
  PyDelphesConfReader() = default;
  explicit PyDelphesConfReader(const pybind11::dict &);

  void ReadFile(std::string_view fileName) override;
  const DelphesParameters &Parameters() const override { return fParams; }

private:
  DelphesParameters ParseDict(const pybind11::dict &) const;

  DelphesParameters fParams;
  //pybind11::scoped_interpreter fInterpreter{};
};

#endif // DelphesPython_PyDelphesConfig_h
