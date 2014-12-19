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
#include <set>

class TObjArray;
class Candidate;

class ExRootTreeBranch;

class DelphesFactory: public TNamed
{
public:
  
  DelphesFactory(const char *name = "ObjectFactory");
  ~DelphesFactory();

  void Clear();
 
  TObjArray *NewPermanentArray();

  TObjArray *NewArray() { return New<TObjArray>(); }

  Candidate *NewCandidate();

  TObject *New(TClass *cl);

  template<typename T>
  T *New() { return static_cast<T *>(New(T::Class())); }

private:

  ExRootTreeBranch *fObjArrays; //!

#if !defined(__CINT__) && !defined(__CLING__)
  std::map< const TClass*, ExRootTreeBranch* > fBranches; //!
#endif

  std::set< TObject* > fPool; //!
  
  ClassDef(DelphesFactory, 1)
};

#endif /* DelphesFactory */

