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
 *  Class handling creation of Candidate,
 *  TObjArray and all other objects.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "TNamed.h"

#include <map>
#include <vector>

class Candidate;

class ExRootTreeBranch;

class DelphesFactory: public TNamed
{
public:
  DelphesFactory(const char *name = "ObjectFactory");

  virtual void Clear(Option_t *option = "");

  bool Has(std::string_view collectionName) const;
  template <typename T>
  std::shared_ptr<T> Book(std::string_view collectionName)
  {
    const std::string name{collectionName};
    if(fMemorySlots.count(name) == 0)
      fMemorySlots[name] = reinterpret_cast<void *>(new T);
    return Attach<T>(collectionName);
  }
  template <typename T>
  std::shared_ptr<T> Attach(std::string_view collectionName)
  {
    if(!Has(collectionName)) throwAttachingFailure(collectionName);
    return std::shared_ptr<T>(reinterpret_cast<T *>(fMemorySlots[std::string{collectionName}]));
  }

  Candidate *NewCandidate();

private:
  void throwAttachingFailure(std::string_view collectionName) const;

  std::map<std::string, void *> fMemorySlots;
  std::vector<std::unique_ptr<TObject> > fCandidates;

  ClassDef(DelphesFactory, 2)
};

#endif /* DelphesFactory */
