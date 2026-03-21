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
  fDelphesFactory(std::make_unique<DelphesFactory>())
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

  const ExRootConfReader::ExRootTaskMap *modules = confReader->GetModules();

  gRandom->SetSeed(confReader->GetInt("::RandomSeed", 0));

  ExRootConfParam param = confReader->GetParam("::ExecutionPath");
  for(int i = 0; i < param.GetSize(); ++i)
  {
    TString name = param[i].GetString(); // retrieve the module name
    if(modules->count(name) == 0)
    {
      std::ostringstream message;
      message << "module '" << name;
      message << "' is specified in ExecutionPath but not configured.";
      throw std::runtime_error(message.str());
    }
    const std::string &delphesModuleName = modules->at(name).Data();
    std::unique_ptr<DelphesModule> moduleObject = DelphesProcessingModuleFactory::Get().Build(delphesModuleName);
    moduleObject->SetFactory(GetFactory());
    moduleObject->SetConfReader(GetConfReader());
    moduleObject->SetName(name);
    if(moduleObject->IsWriter())
    {
      DelphesWriter *writerModule = static_cast<DelphesWriter *>(moduleObject.get());
      writerModule->SetOutputFile(GetOutputFile());
    }
    fModules.emplace_back(std::make_pair(name, std::move(moduleObject)));
  }
}

//------------------------------------------------------------------------------

void Delphes::InitTask()
{
  Init();
  for(const auto &[moduleName, moduleObject] : fModules)
  {
    std::cout << std::left;
    std::cout << std::setw(30) << "** INFO: initializing module";
    std::cout << std::setw(25) << moduleName << std::endl;
    moduleObject->Init();
  }
}

//------------------------------------------------------------------------------

void Delphes::ProcessTask()
{
  for(const auto &[moduleName, moduleObject] : fModules)
    moduleObject->Process();
}

//------------------------------------------------------------------------------

void Delphes::FinishTask()
{
  for(const auto &[moduleName, moduleObject] : fModules)
    moduleObject->Finish();
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
