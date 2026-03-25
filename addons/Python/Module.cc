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
#include <pybind11/stl.h>

#include "classes/DelphesFactory.h"
#include "classes/DelphesReader.h"

#include "PyDelphes.h"
#include "PyDelphesConfReader.h"

#include <dlfcn.h>

#include <ExRootAnalysis/ExRootProgressBar.h>

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

  // on the other hand, Reader are really Python objects
  py::class_<DelphesReader>(m, "Reader", "External event reader object")
    .def(py::init([](std::string_view readerType, const py::kwargs &readerArgs) {
      std::unique_ptr<DelphesReader> readerObj = DelphesReaderFactory::Get().Build(std::string{readerType}, PyDelphesConfReader{readerArgs}.Parameters());
      if(readerArgs.contains("inputFiles"))
      {
        if(const std::vector<std::string> inputFiles = readerArgs["inputFiles"].cast<std::vector<std::string> >(); !inputFiles.empty())
        {
          readerObj->LoadInputFile(inputFiles.at(0));
          readerObj->Clear();
        }
        //TODO: manage multi-files inputs
      }
      return readerObj;
    }));

  py::class_<PyDelphesEvent>(m, "_DelphesEvent", "Event collections content");

  // this is where all the magic is done
  py::class_<PyDelphes>(m, "Delphes", "Main Delphes processing module")
    .def(py::init<>())
    .def_property("reader", &PyDelphes::GetReader, &PyDelphes::SetReader, "Reader module used in event consumption")
    .def_property("modules", &PyDelphes::GetModules, &PyDelphes::SetModules, "Processing modules chain definition")
    .def("next", &PyDelphes::Next, "Perform a new event readout and processing");
}

//------------------------------------------------------------------------------
