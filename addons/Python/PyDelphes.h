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

/** \class PyDelphes' module
 *
 *  Binding module for Python
 *
 *  \author L. Forthomme - AGH, Kraków
 *
 */

#ifndef Python_PyDelphes_h
#define Python_PyDelphes_h

#include "classes/DelphesWriter.h"
#include "modules/Delphes.h"

namespace pybind11
{
class dict;
} // namespace pybind11

class DelphesFactory;

class PyDelphesEvent
{
public:
  explicit PyDelphesEvent(const DelphesFactory &factory) : fFactory(factory) {}
  CandidatesCollection Get() const;

private:
  const DelphesFactory &fFactory;
};

class PyDelphes: public Delphes
{
public:
  PyDelphes() {}
  ~PyDelphes();

  void SetReader(DelphesReader *readerObj) override;

  void Init() override;
  const PyDelphesEvent &Next();
  pybind11::dict GetModules() const;
  void SetModules(const pybind11::dict &modulesObj);

private:
  std::unique_ptr<DelphesReader> fEventReader;
  std::unique_ptr<PyDelphesEvent> fDelphesEvent;
  DelphesParameters fConfig;
  bool fIsInitialised{false};
};

#endif
