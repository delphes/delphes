/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2026  Universite catholique de Louvain (UCL), Belgium
 *                           AGH University of Krakow, Poland
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

/** \class PyDelphes' module
 *
 *  Binding module for Python
 *
 *  \author L. Forthomme - AGH, Kraków
 *
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "classes/DelphesFactory.h"
#include "classes/DelphesReader.h"
#include "classes/DelphesWriter.h"
#include "modules/Delphes.h"

#include "PyDelphes.h"
#include "PyDelphesEvent.h"
#include "PyDelphesParameters.h"

#include <ExRootAnalysis/ExRootProgressBar.h>

namespace py = pybind11;

PyDelphes::~PyDelphes() { FinishTask(); }

//------------------------------------------------------------------------------

void PyDelphes::Init()
{
  ClearModules();
  const DelphesParameters fullConfig = fConfig.cast<DelphesParameters>();
  DelphesFactory *factory = GetFactory();
  if(!factory)
    throw std::runtime_error("Delphes factory was not initialised for this PyDelphes module.");
  fDelphesEvent = std::make_unique<PyDelphesEvent>(*factory);
  for(const std::string &moduleName : fullConfig.Get<std::vector<std::string> >("ExecutionPath"))
  {
    if(!fullConfig.Has<DelphesParameters>(moduleName))
    {
      std::ostringstream message;
      message << "module '" << moduleName;
      message << "' is specified in ExecutionPath but not configured.";
      throw std::runtime_error(message.str());
    }
    try
    {
      const DelphesParameters moduleParams = fullConfig.Get<DelphesParameters>(moduleName);
      const std::string moduleTypeFromParams = moduleParams.Get<std::string>("ModuleType", moduleName);
      std::unique_ptr<DelphesModule> moduleObject = DelphesProcessingModuleFactory::Get().Build(moduleTypeFromParams, moduleParams);
      moduleObject->SetName(moduleName);
      moduleObject->SetFactory(factory);
      if(moduleObject->IsWriter())
      {
        std::cout << "output file: " << GetOutputFile() << std::endl;
        DelphesWriter *writerModule = static_cast<DelphesWriter *>(moduleObject.get());
        writerModule->SetOutputFile(GetOutputFile());
      }
      AddModule(moduleName, moduleObject);
    }
    catch(const std::runtime_error &error)
    {
      std::ostringstream message;
      message << "Failed to build '" << moduleName << "' module. Error: " << error.what();
      if(fullConfig.Has<DelphesParameters>(moduleName))
        message << "\nParameters:\n"
                << fullConfig.Get<DelphesParameters>(moduleName);
      throw std::runtime_error(message.str());
    }
  }
}

//------------------------------------------------------------------------------

void PyDelphes::SetReaderObj(DelphesReader *readerObj)
{
  fEventReader.reset(readerObj);
  //fEventReader = readerObj;
  //Delphes::SetReader(readerObj);
  fEventReader->SetFactory(GetFactory());
  const DelphesParameters fullConfig = fConfig.cast<DelphesParameters>();
  fEventReader->SetMaxEvents(fullConfig.Get<int>("MaxEvents", 0));
  fEventReader->SetSkipEvents(fullConfig.Get<int>("SkipEvents", 0));
}

//------------------------------------------------------------------------------

void PyDelphes::SetReaderConfig(const pybind11::dict &readerArgs)
{
  fReaderConfig = readerArgs;
  const DelphesParameters fullConfig = fReaderConfig.cast<DelphesParameters>();
  fEventReader = DelphesReaderFactory::Get().Build(fullConfig.Get<std::string>("ReaderType"), fullConfig);
  if(readerArgs.contains("inputFiles"))
  {
    if(const std::vector<std::string> inputFiles = readerArgs["inputFiles"].cast<std::vector<std::string> >(); !inputFiles.empty())
      fEventReader->LoadInputFile(inputFiles.at(0));
    //TODO: manage multi-files inputs
  }
  fEventReader->SetFactory(GetFactory());
  fEventReader->SetMaxEvents(fullConfig.Get<int>("MaxEvents", 0));
  fEventReader->SetSkipEvents(fullConfig.Get<int>("SkipEvents", 0));
}

//------------------------------------------------------------------------------

const PyDelphesEvent &PyDelphes::Next()
{
  if(!fEventReader)
    throw std::runtime_error("Event reader object is not yet initialised.");
  if(!fIsInitialised || !fConfig.is(fLastProcessingConfig))
  {
    std::cout << "Configuration has changed!!" << std::endl;
    std::cout << fConfig.cast<DelphesParameters>() << std::endl;
    InitTask();
    fLastProcessingConfig = fConfig;
    fIsInitialised = true;
  }
  Clear();
  fEventReader->ReadEvent();
  ProcessTask();
  Clear();
  fEventReader->Clear();
  return *fDelphesEvent;
}

//------------------------------------------------------------------------------

const py::dict &PyDelphes::GetModules() const
{
  return fConfig;
}

//------------------------------------------------------------------------------

void PyDelphes::SetModules(const py::dict &modulesObj)
{
  py::list execPath;
  for(const std::pair<py::handle, py::handle> &moduleObj : modulesObj)
  {
    const std::string moduleName = moduleObj.first.cast<std::string>();
    if(py::isinstance<py::dict>(moduleObj.second))
    {
      fConfig[py::str(moduleName)] = moduleObj.second.cast<py::dict>();
      execPath.append(py::str(moduleName));
    }
  }
  fConfig[py::str("ExecutionPath")] = execPath;
}

//------------------------------------------------------------------------------
