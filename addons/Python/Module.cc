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

#include <pybind11/pybind11.h>

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesReader.h"
#include "modules/Delphes.h"

#include <dlfcn.h>

namespace py = pybind11;

class PyDelphes: public Delphes
{
public:
  PyDelphes() {}
  ~PyDelphes() { FinishTask(); }

  void Next()
  {
    if(!fIsInitialised)
    {
      InitTask();
      fIsInitialised = true;
    }
    GetReader()->ReadEvent();
    ProcessTask();
    Clear();
    GetReader()->Clear();
  }

private:
  bool fIsInitialised{false};
};

PYBIND11_MODULE(DelphesPython, m)
{
  m.doc() = R"pbdoc(
    Delphes Python plugin
    ---------------------

    A Python wrapper for the processing of events through the Delphes framework.
  )pbdoc";

  m.def("load", [](std::string libraryName) {
    if(::dlopen(libraryName.data(), RTLD_LAZY | RTLD_GLOBAL) == nullptr)
    {
      std::ostringstream message;
      message << "Failed to load library " << libraryName << ".";
      if(const char *err = dlerror(); err != nullptr) message << "\n\t" << err;
      throw std::runtime_error(message.str());
    } }, "Load an external library into the runtime environment");

  py::class_<PyDelphes>(m, "Delphes")
    .def(py::init<>())
    .def_property("reader", &PyDelphes::GetReader, &PyDelphes::SetReader, "Reader module used in event consumption")
    .def("next", &PyDelphes::Next, "Perform a new event readout and processing");
  py::class_<DelphesReader>(m, "Reader")
    .def(py::init([](std::string name, const py::kwargs &) { return DelphesReaderFactory::Get().Build(name); }));
  py::class_<DelphesModule>(m, "Module")
    .def(py::init([](std::string name, const py::kwargs &) { return DelphesProcessingModuleFactory::Get().Build(name); }));
}
