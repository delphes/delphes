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
#include "classes/DelphesReader.h"
#include "classes/DelphesWriter.h"

#include <ExRootAnalysis/ExRootConfReader.h>

#include <TRandom3.h>

Delphes::Delphes(const char *name) :
  fDelphesFactory(std::make_unique<DelphesFactory>("ObjectFactory"))
{
  SetName(name);
}

//------------------------------------------------------------------------------

void Delphes::Clear(Option_t * /*option*/)
{
  if(fDelphesFactory) fDelphesFactory->Clear();
}

//------------------------------------------------------------------------------

void Delphes::Init()
{
  if(!fReader)
    throw std::runtime_error("Failed to initialise the main Delphes module with no reader declared.");

  ExRootConfReader *confReader = GetConfReader();
  confReader->SetName("ConfReader");

  const ExRootConfReader::ExRootTaskMap *modules = confReader->GetModules();
  ExRootConfReader::ExRootTaskMap::const_iterator itModules;

  gRandom->SetSeed(confReader->GetInt("::RandomSeed", 0));

  ExRootConfParam param = confReader->GetParam("::ExecutionPath");
  for(int i = 0; i < param.GetSize(); ++i)
  {
    TString name = param[i].GetString(); // retrieve the module name
    itModules = modules->find(name);
    if(itModules != modules->end()) // found the module in the execution path
    {
      std::unique_ptr<DelphesModule> moduleObject = DelphesProcessingModuleFactory::Get().Build(itModules->second.Data());
      moduleObject->SetFactory(GetFactory());
      if(moduleObject->IsWriter())
      {
        DelphesWriter *writerModule = static_cast<DelphesWriter *>(moduleObject.get());
        writerModule->SetOutputFile(GetOutputFile());
      }
      if(ExRootTask *task = NewTask(dynamic_cast<ExRootTask *>(moduleObject.release()), itModules->first); task)
        Add(task);
    }
    else
    {
      std::ostringstream message;
      message << "module '" << name;
      message << "' is specified in ExecutionPath but not configured.";
      throw std::runtime_error(message.str());
    }
  }
}

//------------------------------------------------------------------------------

void Delphes::SetReader(DelphesReader *reader)
{
  fReader = reader;
  fReader->SetFactory(GetFactory());
}

//------------------------------------------------------------------------------

DelphesFactory *Delphes::GetFactory() const { return fDelphesFactory.get(); }

//------------------------------------------------------------------------------
