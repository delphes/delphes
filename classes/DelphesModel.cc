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

#include "DelphesModel.h"

DelphesModel::DelphesModel(std::string_view model_name) : name_(model_name), model_(RNTupleModel::Create().release()) {}

DelphesModel::~DelphesModel() {}

const REntry &DelphesModel::Entry() const
{
  if(!entry_)
  {
    std::cout << "Entry was not specified in module, using the default one from '" << name_ << "' model." << std::endl;
    if(!model_)
      throw std::runtime_error(std::string{"Base model for '"} + name_ + "' payloads collection was already released.");
    return model_->GetDefaultEntry();
  }
  return *entry_;
}

REntry &DelphesModel::Entry()
{
  if(!entry_)
  {
    std::cout << "Non-const-qualified '" << name_ << "' entry is nonexistent. Building a new one from the model." << std::endl;
    if(!model_)
      throw std::runtime_error(std::string{"Base model for '"} + name_ + "' payloads collection was already released.");
    entry_ = model_->CreateEntry().release();
  }
  return *entry_;
}

void DelphesModel::throwBookingFailure(std::string_view collection_name, std::string_view details) const noexcept(false)
{
  std::ostringstream os;
  os << "Failed to book '" << name_ << "' collection with name '" << collection_name
     << "'." + (!details.empty() ? std::string{"\nDetails: "} + std::string{details} : std::string{});
  throw std::runtime_error(os.str());
}

void DelphesModel::throwAttachingFailure(std::string_view collection_name, std::string_view details) const noexcept(false)
{
  std::stringstream os;
  os << "Failed to attach memory segment to '" << name_ << "' collection with '"
     << collection_name << "' label.\n"
     << "List of fields registered:";
  for(const auto &value : Entry())
    os << "\n - '" << value.GetField().GetFieldName() << "' of type '" << value.GetField().GetTypeName() << "'";
  if(!details.empty())
    os << "\nDetails: " << details;
  throw std::runtime_error(os.str());
}
