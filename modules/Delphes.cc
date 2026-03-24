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
#include "classes/DelphesTCLConfReader.h"
#include "classes/DelphesWriter.h"

#include <ExRootAnalysis/ExRootProgressBar.h>

#include <TRandom3.h>

Delphes::Delphes(const char *name) :
  DelphesModule(DelphesParameters{}),
  fDelphesFactory(std::make_unique<DelphesFactory>())
{
  SetName(name);
}

//------------------------------------------------------------------------------

void Delphes::Clear()
{
  if(fDelphesFactory) fDelphesFactory->Clear();
}

//------------------------------------------------------------------------------

void Delphes::Init()
{
  if(!fReader)
    throw std::runtime_error("Failed to initialise the main Delphes module with no reader declared.");
  if(!fConfReader)
    throw std::runtime_error("Failed to initialise the main Delphes module with no user configuration reader declared.");

  const auto userConfig = fConfReader->Parameters();
  gRandom->SetSeed(userConfig.Get<int>("RandomSeed", 0));

  for(const std::string &moduleName : userConfig.Get<std::vector<std::string> >("ExecutionPath"))
  {
    if(!userConfig.Has<DelphesParameters>(moduleName))
    {
      std::ostringstream message;
      message << "module '" << moduleName;
      message << "' is specified in ExecutionPath but not configured.";
      throw std::runtime_error(message.str());
    }
    try
    {
      const DelphesParameters moduleParams = userConfig.Get<DelphesParameters>(moduleName);
      const std::string moduleTypeFromParams = moduleParams.Get<std::string>("ModuleType", moduleName);
      std::unique_ptr<DelphesModule> moduleObject = DelphesProcessingModuleFactory::Get().Build(moduleTypeFromParams, moduleParams);
      moduleObject->SetName(moduleName);
      moduleObject->SetFactory(GetFactory());
      if(moduleObject->IsWriter())
      {
        DelphesWriter *writerModule = static_cast<DelphesWriter *>(moduleObject.get());
        writerModule->SetOutputFile(GetOutputFile());
      }
      fModules.emplace_back(std::make_pair(moduleName, std::move(moduleObject)));
    }
    catch(const std::runtime_error &error)
    {
      std::ostringstream message;
      message << "Failed to build '" << moduleName << "' module. Error: " << error.what();
      if(userConfig.Has<DelphesParameters>(moduleName))
        message << "\nParameters:\n"
                << userConfig.Get<DelphesParameters>(moduleName);
      throw std::runtime_error(message.str());
    }
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
