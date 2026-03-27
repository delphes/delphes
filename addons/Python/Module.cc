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

#include <pybind11/stl.h>

#include "classes/DelphesFactory.h"
#include "classes/DelphesModule.h"
#include "classes/DelphesReader.h"

#include "PyDelphes.h"
#include "PyDelphesEvent.h"

#include <dlfcn.h>

#include <ExRootAnalysis/ExRootProgressBar.h>
#include <TStopwatch.h>

namespace py = pybind11;

//------------------------------------------------------------------------------

PYBIND11_MODULE(DelphesPython, m)
{
  m.doc() = R"pbdoc(
    Delphes Python plugin
    ---------------------

    A Python wrapper for the processing of events through the Delphes framework.
  )pbdoc";

  m.def("load", [](std::string_view libraryName) {
    if(::dlopen(libraryName.data(), RTLD_LAZY | RTLD_GLOBAL) == nullptr)
    {
      std::ostringstream message;
      message << "Failed to load library '" << libraryName << "'.";
      if(const char *err = dlerror(); err != nullptr) message << " " << err;
      throw py::import_error(message.str());
    }
    py::print("Successfully loaded the library", libraryName, "into the runtime environment."); }, "Load an external library into the runtime environment");

  // mimic Module as an actual Python object, while it is only a translator for a Python dictionary
  m.def("Module", [](std::string_view moduleType, const py::kwargs &moduleArgs) -> py::dict {
    py::dict paramsObj = moduleArgs;
    paramsObj[py::str("ModuleType")] = py::str(moduleType);
    return paramsObj; }, "Define a processing module configuration");

  // mimic Reader as an actual Python object, while it is only a translator for a Python dictionary
  m.def("Reader", [](std::string_view moduleType, const py::kwargs &readerArgs) -> py::dict {
    py::dict paramsObj = readerArgs;
    paramsObj[py::str("ReaderType")] = py::str(moduleType);
    return paramsObj; }, "Define an external event reader configuration");

  py::class_<PyDelphesEvent>(m, "_DelphesEvent", py::dynamic_attr(), "Event collections content")
    .def("__getitem__", &PyDelphesEvent::Get<std::vector<Candidate *> >, "Retrieve a collection from the event");

  // this is where all the magic is done
  py::class_<PyDelphes>(m, "Delphes", "Main Delphes processing module")
    .def(py::init<>())
    .def_property("reader", &PyDelphes::GetReaderConfig, &PyDelphes::SetReaderConfig, "Reader module used in event consumption")
    .def_property("modules", &PyDelphes::GetModules, &PyDelphes::SetModules, "Processing modules chain definition")
    .def("next", &PyDelphes::Next, "Perform a new event readout and processing");

  py::class_<Candidate, std::unique_ptr<Candidate, py::nodelete> >(m, "Candidate")
    .def(py::init<>())
    .def_property_readonly("pt", [](Candidate *self) {
      std::cout << "--> " << self << std::endl;
      if (!self) return -1.f;
      return self->PT; });
  //.def_readwrite("pt", &Candidate::PT);
}

//------------------------------------------------------------------------------
