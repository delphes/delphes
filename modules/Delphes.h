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

#ifndef Delphes_h
#define Delphes_h

/** \class Delphes
 *
 *  Main Delphes module.
 *  Controls execution of all other modules.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TFolder;
class TObjArray;

class ExRootTreeWriter;

class DelphesFactory;

class Delphes: public DelphesModule
{
public:

  Delphes(const char *name = "Delphes");
  ~Delphes();

  void SetTreeWriter(ExRootTreeWriter *treeWriter);
  
  DelphesFactory *GetFactory() const { return fFactory; }

  void Clear();

  virtual void Init();
  virtual void Process();
  virtual void Finish();

private:

  DelphesFactory *fFactory;

  ClassDef(Delphes, 1)
};

#endif /* Delphes_h */

