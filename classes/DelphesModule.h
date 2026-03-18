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

#ifndef DelphesModule_h
#define DelphesModule_h

/** \class DelphesModule
 *
 *  Base class for all Delphes modules
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootTask.h"

class TClass;
class TObject;
class TFolder;
class TClonesArray;

class Candidate;

class ExRootResult;
class ExRootTreeBranch;
class ExRootTreeWriter;

class DelphesFactory;

class DelphesModule: public ExRootTask
{
public:
  DelphesModule();
  virtual ~DelphesModule();

  virtual void Init();
  virtual void Process();
  virtual void Finish();

  std::shared_ptr<std::vector<Candidate *> > &ImportArray(const char *name);
  std::shared_ptr<std::vector<Candidate *> > &ExportArray(const char *name);

  ExRootTreeBranch *NewBranch(const char *name, TClass *cl);
  void AddInfo(const char *name, Double_t value);

  ExRootResult *GetPlots();
  DelphesFactory *GetFactory();

protected:
  ExRootTreeWriter *fTreeWriter{nullptr};
  DelphesFactory *fFactory{nullptr};

private:
  ExRootResult *fPlots{nullptr};

  TFolder *fPlotFolder{nullptr}, *fExportFolder{nullptr};

  ClassDef(DelphesModule, 1)
};

#endif /* DelphesModule_h */
