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

#ifndef DelphesFactory_h
#define DelphesFactory_h

/** \class DelphesFactory
 *
 *  Class handling creation of Candidate, and all other objects.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *  \author L. Forthomme - AGH, Krakow
 *
 */

#include <map>
#include <memory>
#include <string>
#include <vector>

class Candidate;

class DelphesFactory
{
public:
  DelphesFactory() = default;

  virtual void Clear(); ///< Clear all collections booked in this factory

  bool Has(std::string_view collectionName) const; ///< Is the collection already booked?
  /// Book the memory segment for a new collection
  template <typename T>
  std::shared_ptr<T> Book(std::string_view collectionName, bool forExport)
  {
    const std::string name{collectionName};
    if(fMemorySlots.count(name) == 0)
      fMemorySlots[name] = reinterpret_cast<void *>(new T);
    if(forExport)
      fExportCollections.emplace_back(std::make_pair(name, &typeid(T)));
    return Attach<T>(collectionName);
  }
  /// Attach a pointer to a collection handled by this factory
  template <typename T>
  std::shared_ptr<T> Attach(std::string_view collectionName) const
  {
    if(!Has(collectionName)) ThrowAttachingFailure(collectionName);
    return std::shared_ptr<T>(reinterpret_cast<T *>(fMemorySlots.at(std::string{collectionName})), [](T *) {});
  }
  void *Attach(std::string_view collectionName) const
  {
    if(!Has(collectionName)) ThrowAttachingFailure(collectionName);
    return fMemorySlots.at(std::string{collectionName});
  }
  std::vector<std::string> GetCollections() const; ///< Retrieve the name of all collections booked
  /// Retrieve the list of collections to be exported automatically
  const std::vector<std::pair<std::string, const std::type_info *> > &GetExportCollections() const { return fExportCollections; }

  Candidate *NewCandidate(); ///< Construct a new candidate to fill a collection

private:
  void ThrowAttachingFailure(std::string_view collectionName) const;

  std::map<std::string, void *> fMemorySlots;
  std::vector<std::pair<std::string, const std::type_info *> > fExportCollections;
  std::vector<std::unique_ptr<Candidate> > fCandidates;
};

#endif /* DelphesFactory */
