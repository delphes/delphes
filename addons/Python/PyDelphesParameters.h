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

#ifndef DelphesPython_PyDelphesParameters_h
#define DelphesPython_PyDelphesParameters_h

#include "classes/DelphesParameters.h"

#include <pybind11/pybind11.h>

class PyDelphes;

namespace pybind11
{
class dict;
} // namespace pybind11

namespace pybind11
{
namespace detail
{

template <>
struct type_caster<DelphesParameters>
{
  PYBIND11_TYPE_CASTER(DelphesParameters, _("dict"));

  static handle
  cast(const DelphesParameters &paramObj, return_value_policy /*policy*/, handle /*parent*/)
  {
    pybind11::dict pyDict;
    for(const std::string &paramKey : paramObj.Keys())
    {
      pybind11::str pyKey(paramKey);
      if(paramObj.Has<DelphesParameters>(paramKey))
        pyDict[pyKey] = pybind11::cast(paramObj.Get<DelphesParameters>(paramKey));
      else if(paramObj.Has<int>(paramKey))
        pyDict[pyKey] = pybind11::cast(paramObj.Get<int>(paramKey));
      else if(paramObj.Has<double>(paramKey))
        pyDict[pyKey] = pybind11::cast(paramObj.Get<double>(paramKey));
      else if(paramObj.Has<std::vector<int> >(paramKey))
        pyDict[pyKey] = pybind11::cast(paramObj.Get<std::vector<int> >(paramKey));
      else if(paramObj.Has<std::vector<DelphesParameters> >(paramKey))
      {
        pybind11::list pyList;
        for(const DelphesParameters &subParams : paramObj.Get<std::vector<DelphesParameters> >(paramKey))
          pyList.append(pybind11::cast(subParams));
        pyDict[pyKey] = pyList;
      }
      else if(paramObj.Has<std::vector<double> >(paramKey))
        pyDict[pyKey] = pybind11::cast(paramObj.Get<std::vector<double> >(paramKey));
      else if(paramObj.Has<std::vector<std::vector<std::string> > >(paramKey))
        pyDict[pyKey] = pybind11::cast(paramObj.Get<std::vector<std::vector<std::string> > >(paramKey));
      else if(paramObj.Has<std::vector<std::vector<DelphesParameters> > >(paramKey))
      {
        pybind11::list pyList;
        for(const std::vector<DelphesParameters> &subParams :
          paramObj.Get<std::vector<std::vector<DelphesParameters> > >(paramKey))
        {
          pybind11::list subPyList;
          for(const DelphesParameters &subSubParams : subParams)
            subPyList.append(pybind11::cast(subSubParams));
          pyList.append(subPyList);
        }
        pyDict[pyKey] = pyList;
      }
      else if(paramObj.Has<std::vector<std::vector<int> > >(paramKey))
        pyDict[pyKey] = pybind11::cast(paramObj.Get<std::vector<std::vector<int> > >(paramKey));
      else if(paramObj.Has<std::vector<std::vector<double> > >(paramKey))
        pyDict[pyKey] = pybind11::cast(paramObj.Get<std::vector<std::vector<double> > >(paramKey));
      else if(paramObj.Has<std::vector<std::string> >(paramKey))
        pyDict[pyKey] = pybind11::cast(paramObj.Get<std::vector<std::string> >(paramKey));
      else if(paramObj.Has<std::string>(paramKey))
        pyDict[pyKey] = pybind11::cast(paramObj.Get<std::string>(paramKey));
      /*else
        std::cout << "** WARNING: failed to parse the Delphes parameters with key '" << paramKey << "' and YAML representation:\n"
                  << paramObj.Get<std::string>(paramKey);*/
    }
    return pyDict.release();
  }

  bool load(handle src, bool /*convert*/)
  {
    if(!pybind11::isinstance<pybind11::dict>(src))
      return false;
    pybind11::dict dictObj = pybind11::reinterpret_borrow<pybind11::dict>(src);
    for(const std::pair<pybind11::handle, pybind11::handle> &dictItem : dictObj)
    {
      const std::string itemKey = dictItem.first.cast<std::string>();
      if(pybind11::isinstance<pybind11::dict>(dictItem.second))
        value.Set(itemKey, dictItem.second.cast<DelphesParameters>());
      else if(pybind11::isinstance<pybind11::list>(dictItem.second))
      {
        const pybind11::handle &firstItem = dictItem.second.cast<pybind11::list>()[0];
        if(pybind11::isinstance<pybind11::dict>(firstItem))
        {
          std::vector<DelphesParameters> paramsList;
          for(const pybind11::handle &listItem : dictItem.second.cast<pybind11::list>())
            paramsList.emplace_back(listItem.cast<DelphesParameters>());
          value.Set(itemKey, paramsList);
        }
        else if(pybind11::isinstance<pybind11::int_>(firstItem))
          value.Set(itemKey, dictItem.second.cast<std::vector<int> >());
        else if(pybind11::isinstance<pybind11::float_>(firstItem))
          value.Set(itemKey, dictItem.second.cast<std::vector<double> >());
        else if(pybind11::isinstance<pybind11::str>(firstItem))
          value.Set(itemKey, dictItem.second.cast<std::vector<std::string> >());
        else if(pybind11::isinstance<pybind11::list>(firstItem) && !firstItem.cast<pybind11::list>().empty()) // list of lists
        {
          const pybind11::handle &firstFirstItem = firstItem.cast<pybind11::list>()[0];
          if(pybind11::isinstance<pybind11::int_>(firstFirstItem))
            value.Set(itemKey, dictItem.second.cast<std::vector<std::vector<int> > >());
          else if(pybind11::isinstance<pybind11::float_>(firstFirstItem))
            value.Set(itemKey, dictItem.second.cast<std::vector<std::vector<double> > >());
          else if(pybind11::isinstance<pybind11::str>(firstFirstItem))
            value.Set(itemKey, dictItem.second.cast<std::vector<std::vector<std::string> > >());
          else
            value.Set(itemKey, std::vector<std::vector<std::string> >{});
        }
      }
      else if(pybind11::isinstance<pybind11::int_>(dictItem.second))
        value.Set(itemKey, dictItem.second.cast<int>());
      else if(pybind11::isinstance<pybind11::float_>(dictItem.second))
        value.Set(itemKey, dictItem.second.cast<double>());
      else if(pybind11::isinstance<pybind11::str>(dictItem.second))
        value.Set(itemKey, dictItem.second.cast<std::string>());
      else
        return false;
    }
    return true;
  }
};

} // namespace detail
} // namespace pybind11

class PyDelphesParameters
{
public:
  explicit PyDelphesParameters(PyDelphes &moduleObj);

  void set_item(std::string_view key, pybind11::handle value);
  pybind11::object get_item(std::string_view key);
  pybind11::object del_item(std::string_view key);

  size_t len();

  pybind11::iterator iter();

  bool contains(std::string_view key);

  pybind11::list keys();
  pybind11::list values();
  pybind11::list items();

  void update(const pybind11::dict &other);
  void clear();

  pybind11::object pop(std::string_view key, pybind11::object default_val = pybind11::none());

private:
  void sync();

  PyDelphes &fModule;
  pybind11::dict fParams;
};

#endif // DelphesPython_PyDelphesParameters_h
