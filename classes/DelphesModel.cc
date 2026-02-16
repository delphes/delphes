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

DelphesModel::DelphesModel(std::string_view model_name) :
  fName(model_name), fModel(RNTupleModel::Create().release()) {}

DelphesModel::~DelphesModel() {}

const REntry &DelphesModel::Entry() const
{
  if(!fEntry)
  {
    std::cout << "Entry was not specified in module, using the default one from '" << fName << "' model." << std::endl;
    if(!fModel)
      throw std::runtime_error(std::string{"Base model for '"} + fName + "' fields collection was already released.");
    return fModel->GetDefaultEntry();
  }
  return *fEntry;
}

REntry &DelphesModel::Entry()
{
  if(!fEntry)
  {
    std::cout << "Non-const-qualified '" << fName << "' entry is nonexistent. Building a new one from the model." << std::endl;
    if(!fModel)
      throw std::runtime_error(std::string{"Base model for '"} + fName + "' fields collection was already released.");
    fEntry = fModel->CreateEntry().release();
  }
  return *fEntry;
}

void DelphesModel::throwBookingFailure(std::string_view field_name, std::string_view details) const noexcept(false)
{
  std::ostringstream os;
  os << "Failed to book '" << fName << "' field with name '" << field_name << "'.";
  if(!details.empty())
    os << "\nDetails: " << details;
  throw std::runtime_error(os.str());
}

void DelphesModel::throwAttachingFailure(std::string_view field_name, std::string_view details) const noexcept(false)
{
  std::stringstream os;
  os << "Failed to attach memory segment to '" << fName << "' field with '" << field_name << "' label.\n"
     << "List of fields registered:";
  for(const auto &value : Entry())
    os << "\n - '" << value.GetField().GetFieldName() << "' of type '" << value.GetField().GetTypeName() << "'";
  if(!details.empty())
    os << "\nDetails: " << details;
  throw std::runtime_error(os.str());
}

DelphesModel::FieldName::FieldName(std::string_view field_label) :
  fLabel(field_label)
{
  // one-liner with C++23
  //fFieldLabel = fLabel | std::views::split("/") | std::views::join_with("__") | std::ranges::to<std::string>();
  fFieldLabel = fLabel;
  auto &&pos = fFieldLabel.find("/", size_t{});
  while(pos != std::string::npos)
  {
    fFieldLabel.replace(pos, 1, "__");
    pos = fFieldLabel.find("/", pos + 2);
  }
}
