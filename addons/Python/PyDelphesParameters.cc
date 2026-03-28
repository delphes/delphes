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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "classes/DelphesClasses.h"
#include "classes/DelphesReader.h"

#include "PyDelphes.h"
#include "PyDelphesEvent.h"
#include "PyDelphesParameters.h"

#include <ExRootAnalysis/ExRootProgressBar.h>

#include <iostream> //FIXME

namespace py = pybind11;

PyDelphesParameters::PyDelphesParameters(PyDelphes &moduleObj) :
  fModule(moduleObj), fParams(moduleObj.GetModulesConfig()) { sync(); }

//------------------------------------------------------------------------------

void PyDelphesParameters::sync()
{
  fModule.SetModulesConfig(fParams);
}

//------------------------------------------------------------------------------

void PyDelphesParameters::set_item(std::string_view key, py::handle value)
{
  fParams[key.data()] = value;
  sync();
}

//------------------------------------------------------------------------------

py::object PyDelphesParameters::get_item(std::string_view key)
{
  return fParams[key.data()];
}

//------------------------------------------------------------------------------

py::object PyDelphesParameters::del_item(std::string_view key)
{
  return fParams.attr("__delitem__")(key);
  sync();
}

//------------------------------------------------------------------------------

size_t PyDelphesParameters::len() { return py::len(fParams); }

//------------------------------------------------------------------------------

py::iterator PyDelphesParameters::iter() { return fParams.attr("__iter__")(); }

//------------------------------------------------------------------------------

bool PyDelphesParameters::contains(std::string_view key)
{
  return fParams.contains(py::str(key));
}

//------------------------------------------------------------------------------

py::list PyDelphesParameters::keys() { return fParams.attr("keys")(); }

//------------------------------------------------------------------------------

py::list PyDelphesParameters::values() { return fParams.attr("values")(); }

//------------------------------------------------------------------------------

py::list PyDelphesParameters::items() { return fParams.attr("items")(); }

//------------------------------------------------------------------------------

void PyDelphesParameters::update(const py::dict &other)
{
  fParams.attr("update")(other);
  sync();
}

//------------------------------------------------------------------------------

void PyDelphesParameters::clear()
{
  fParams.clear();
  sync();
}

//------------------------------------------------------------------------------

py::object PyDelphesParameters::pop(std::string_view key, py::object default_val)
{
  py::object result = default_val.is_none() ?
    fParams.attr("pop")(key) :
    fParams.attr("pop")(key, default_val);
  sync();
  return result;
}

//------------------------------------------------------------------------------
