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

#include "classes/DelphesModule.h"
#include "classes/DelphesReader.h"

#include "PyDelphes.h"
#include "PyDelphesEvent.h"
#include "PyDelphesParameters.h"

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

  m.def("load", &LoadLibrary, "Load an external library into the runtime environment");

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

  py::class_<Candidate, std::unique_ptr<Candidate, py::nodelete> >(m, "Candidate")
    .def(py::init<>())
    .def_property_readonly("pt", [](const Candidate &self) { return self.Momentum.Pt(); }, "Candidate transverse momentum (in GeV/c)")
    .def_property_readonly("eta", [](const Candidate &self) { return self.Momentum.Eta(); }, "Candidate pseudo-rapidity")
    .def_property_readonly("phi", [](const Candidate &self) { return self.Momentum.Phi(); }, "Candidate azimuthal angle (in rad)")
    .def_property_readonly("vertex", [](const Candidate &self) { return std::array{self.Position.X(), self.Position.Y(), self.Position.Z()}; }, "Candidate vertex position (in m)");

  py::class_<PyDelphesEvent>(m, "_DelphesEvent", py::dynamic_attr(), "Event collections content")
    .def("__getitem__", &PyDelphesEvent::Get<std::vector<Candidate *> >, py::return_value_policy::reference, "Retrieve a collection from the event");

  py::class_<PyDelphesParameters>(m, "PyDelphesParameters")
    .def("__getitem__", &PyDelphesParameters::get_item)
    .def("__setitem__", &PyDelphesParameters::set_item)
    .def("__delitem__", &PyDelphesParameters::del_item)
    .def("__len__", &PyDelphesParameters::len)
    .def("__iter__", &PyDelphesParameters::iter)
    .def("__contains__", &PyDelphesParameters::contains)
    .def("keys", &PyDelphesParameters::keys)
    .def("values", &PyDelphesParameters::values)
    .def("items", &PyDelphesParameters::items)
    .def("update", &PyDelphesParameters::update)
    .def("clear", &PyDelphesParameters::clear)
    .def("pop", &PyDelphesParameters::pop);

  // this is where all the magic is done
  py::class_<PyDelphes>(m, "Delphes", "Main Delphes processing module")
    .def(py::init<>())
    .def_property("reader", &PyDelphes::GetReaderConfig, &PyDelphes::SetReaderConfig, "Reader module used in event consumption")
    .def_property_readonly("modules", [](PyDelphes &self) { return std::make_unique<PyDelphesParameters>(self); }, "Processing modules chain definition")
    .def("loadTCL", &PyDelphes::LoadTCL, "Load configuration from an external TCL file")
    .def("reset", &PyDelphes::Reset, "Reset the reader to its first event")
    .def("next", &PyDelphes::Next, "Perform a new event readout and processing");
}

//------------------------------------------------------------------------------
