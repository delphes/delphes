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

#include <pybind11/stl.h>

#include "PyDelphesConfReader.h"

#include <iostream> //FIXME

namespace py = pybind11;

PyDelphesConfReader::PyDelphesConfReader(const py::dict &dict) { fParams = ParseDict(dict); }

//------------------------------------------------------------------------------

void PyDelphesConfReader::ReadFile(std::string_view fileName)
{
}

//------------------------------------------------------------------------------

DelphesParameters PyDelphesConfReader::ParseDict(const py::dict &dictObj) const
{
  DelphesParameters paramsObj;
  for(const std::pair<py::handle, py::handle> &dictItem : dictObj)
  {
    const std::string itemKey = dictItem.first.cast<std::string>();
    if(py::isinstance<py::dict>(dictItem.second))
      paramsObj.Set(itemKey, ParseDict(dictItem.second.cast<py::dict>()));
    else if(py::isinstance<py::list>(dictItem.second))
    {
      const py::handle &firstItem = dictItem.second.cast<py::list>()[0];
      if(py::isinstance<py::dict>(firstItem))
      {
        std::vector<DelphesParameters> paramsList;
        for(const py::handle &listItem : dictItem.second.cast<py::list>())
          paramsList.emplace_back(ParseDict(listItem.cast<py::dict>()));
        paramsObj.Set(itemKey, paramsList);
      }
      else if(py::isinstance<py::int_>(firstItem))
        paramsObj.Set(itemKey, dictItem.second.cast<std::vector<int> >());
      else if(py::isinstance<py::float_>(firstItem))
        paramsObj.Set(itemKey, dictItem.second.cast<std::vector<double> >());
      else if(py::isinstance<py::str>(firstItem))
        paramsObj.Set(itemKey, dictItem.second.cast<std::vector<std::string> >());
    }
    else if(py::isinstance<py::int_>(dictItem.second))
      paramsObj.Set(itemKey, dictItem.second.cast<int>());
    else if(py::isinstance<py::float_>(dictItem.second))
      paramsObj.Set(itemKey, dictItem.second.cast<double>());
    else if(py::isinstance<py::str>(dictItem.second))
      paramsObj.Set(itemKey, dictItem.second.cast<std::string>());
    else
    {
      std::ostringstream message;
      message << "Failed to unpack the value of parameter '" << itemKey << "'.";
      throw std::runtime_error(message.str());
    }
  }
  return paramsObj;
}

//------------------------------------------------------------------------------
