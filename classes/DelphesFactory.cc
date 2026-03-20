/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
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

/** \class DelphesFactory
 *
 *  Class handling creation of Candidate, and all other objects.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

DelphesFactory::DelphesFactory(const char *name) : TNamed(name, "") {}

//------------------------------------------------------------------------------

void DelphesFactory::Clear(Option_t * /*option*/)
{
  TProcessID::SetObjectCount(0);
  fCandidates.clear();
}

//------------------------------------------------------------------------------

Candidate *DelphesFactory::NewCandidate()
{
  Candidate *object = static_cast<Candidate *>(fCandidates.emplace_back(std::make_unique<Candidate>()).get());
  object->SetFactory(this);
  TProcessID::AssignID(object);
  return object;
}

//------------------------------------------------------------------------------

bool DelphesFactory::Has(std::string_view collectionName) const
{
  return fMemorySlots.count(std::string{collectionName}) > 0;
}

//------------------------------------------------------------------------------

std::vector<std::string> DelphesFactory::GetCollections() const
{
  std::vector<std::string> collections;
  for(const auto &[collectionName, collectionPtr] : fMemorySlots)
    collections.emplace_back(collectionName);
  return collections;
}

//------------------------------------------------------------------------------

void DelphesFactory::throwAttachingFailure(std::string_view collectionName) const
{
  std::ostringstream os;
  os << "Failed to attach memory segment to the collection name '" << collectionName << "'.";
  throw std::runtime_error(os.str());
}

//------------------------------------------------------------------------------
