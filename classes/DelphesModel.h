/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *                2026       AGH University of Krakow, Poland
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

#ifndef DelphesModel_h
#define DelphesModel_h

#include <ROOT/RNTupleModel.hxx>

/** \class DelphesModel
 *
 *  Collection of collections model
 *
 *  \author L. Forthomme - AGH, Krakow
 *  \note Shamelessly (self-)stolen from the TDAnalyser memory management layout
 *
 */

using namespace ROOT; // from ROOT v6.36 on, default namespace
using namespace ROOT::Experimental;
/// Handle to retrieve collections from a data model
template <typename T>
using InputHandle = std::shared_ptr<T>;
/// Handle to store collections into a data model
template <typename T>
using OutputHandle = std::shared_ptr<T>;

class DelphesModel
{
public:
  explicit DelphesModel(std::string_view model_name);
  ~DelphesModel();

  /// Book a memory segment, and attach it to a given field identified by a collection name and a description
  template <typename T>
  OutputHandle<T> book(std::string_view collection_name, std::string_view description = "")
  {
    if(const auto collection_name_str = std::string{collection_name}; !fields_.insert(collection_name_str).second)
      throwBookingFailure(collection_name, "Collection name already exists");
    else if(field_types_.count(collection_name_str) == 0)
      field_types_[collection_name_str] = RField<T>::TypeName();
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 36, 0)
    return OutputHandle<T>{model_->MakeField<T>(collection_name, description)};
#else
    return OutputHandle<T>{model_->MakeField<T>({collection_name, description})};
  };
#endif
  }

  /// Attach a memory segment to a given field identified by a collection name
  template <typename T>
  InputHandle<T> Attach(std::string_view collection_name) const
  {
    try
    {
      return Entry().GetPtr<T>(collection_name);
    }
    catch(const RException &except)
    {
      throwAttachingFailure(collection_name, except.GetError().GetReport());
      throw;
    }
  }

  const REntry &Entry() const;
  REntry &Entry();

private:
  void throwBookingFailure(std::string_view, std::string_view details = "") const noexcept(false);
  void throwAttachingFailure(std::string_view, std::string_view details = "") const noexcept(false);

  const std::string name_; ///< payload collections name (run, event/trigger)
  RNTupleModel *model_{nullptr}; ///< RNTuple data model
  bool model_released_{false}; ///< has the model property been given to a writer?
  REntry *entry_{nullptr}; ///< default entry associated to the data model
  std::set<std::string> fields_; ///< collection of field names managed by this data model
  /// mapping of C++ types associated to each field registered in the data model
  std::unordered_map<std::string, std::string> field_types_;
};

#endif /* DelphesModel_h */
