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

/** \class Delphes
 *
 *  Main Delphes module.
 *  Controls execution of all other modules.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Delphes.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootConfReader.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "TDatabasePDG.h"
#include "TFolder.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <stdio.h>
#include <string.h>

using namespace std;

Delphes::Delphes(const char *name) :
  fFactory(0)
{
  TFolder *folder = new TFolder(name, "");
  fFactory = new DelphesFactory("ObjectFactory");

  SetName(name);
  SetFolder(folder);

  folder->Add(this);
  folder->Add(fFactory);

  gROOT->GetListOfBrowsables()->Add(folder);
}

//------------------------------------------------------------------------------

Delphes::~Delphes()
{
  TFolder *folder = GetFolder();
  if(folder)
  {
    folder->Clear();
    delete folder;
  }
  if(fFactory) delete fFactory;
}

//------------------------------------------------------------------------------

void Delphes::Clear()
{
  if(fFactory) fFactory->Clear();
}

//------------------------------------------------------------------------------

void Delphes::SetTreeWriter(ExRootTreeWriter *treeWriter)
{
  treeWriter->SetName("TreeWriter");
  GetFolder()->Add(treeWriter);
}

//------------------------------------------------------------------------------

void Delphes::Init()
{
  stringstream message;
  ExRootConfReader *confReader = GetConfReader();
  confReader->SetName("ConfReader");
  GetFolder()->Add(confReader);

  TString name;
  ExRootTask *task;
  const ExRootConfReader::ExRootTaskMap *modules = confReader->GetModules();
  ExRootConfReader::ExRootTaskMap::const_iterator itModules;

  ExRootConfParam param = confReader->GetParam("::ExecutionPath");
  Long_t i, size = param.GetSize();

  gRandom->SetSeed(confReader->GetInt("::RandomSeed", 0));

  for(i = 0; i < size; ++i)
  {
    name = param[i].GetString();
    itModules = modules->find(name);
    if(itModules != modules->end())
    {
      task = NewTask(itModules->second, itModules->first);
      if(task)
      {
        task->SetFolder(GetFolder());
        Add(task);
      }
    }
    else
    {
      message << "module '" << name;
      message << "' is specified in ExecutionPath but not configured.";
      throw runtime_error(message.str());
    }
  }
}

//------------------------------------------------------------------------------

void Delphes::Process()
{
}

//------------------------------------------------------------------------------

void Delphes::Finish()
{
}

//------------------------------------------------------------------------------
