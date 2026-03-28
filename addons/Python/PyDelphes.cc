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
#include "classes/DelphesTCLConfReader.h"
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
  const DelphesParameters modulesConfig = fModulesConfig.cast<DelphesParameters>();
  DelphesFactory *factory = GetFactory();
  if(!factory)
    throw std::runtime_error("Delphes factory was not initialised for this PyDelphes module.");
  fDelphesEvent = std::make_unique<PyDelphesEvent>(*factory);
  const std::string outputFile = GetOutputFile();

  if(!fEventReader)
    throw std::runtime_error("Event reader module was not yet initialised at PyDelphes module initialisation.");
  fEventReader->SetMaxEvents(modulesConfig.Get<int>("MaxEvents", 0));
  fEventReader->SetSkipEvents(modulesConfig.Get<int>("SkipEvents", 0));
  fEventReader->Reset();

  for(const std::string &moduleName : modulesConfig.Get<std::vector<std::string> >("ExecutionPath"))
  {
    if(!modulesConfig.Has<DelphesParameters>(moduleName))
    {
      std::ostringstream message;
      message << "module '" << moduleName << "' is specified in ExecutionPath but not configured.";
      throw std::runtime_error(message.str());
    }
    try
    {
      const DelphesParameters moduleParams = modulesConfig.Get<DelphesParameters>(moduleName);
      const std::string moduleTypeFromParams = moduleParams.Get<std::string>("ModuleType", moduleName);
      std::unique_ptr<DelphesModule> moduleObject = DelphesProcessingModuleFactory::Get().Build(moduleTypeFromParams, moduleParams);
      moduleObject->SetName(moduleName);
      moduleObject->SetFactory(factory);
      if(moduleObject->IsWriter())
      {
        if(outputFile.empty()) continue; // skip writers if not output is defined
        std::cout << "output file: " << outputFile << std::endl;
        DelphesWriter *writerModule = static_cast<DelphesWriter *>(moduleObject.get());
        writerModule->SetOutputFile(outputFile);
      }
      AddModule(moduleName, moduleObject);
    }
    catch(const std::runtime_error &error)
    {
      std::ostringstream message;
      message << "Failed to build '" << moduleName << "' module. Error: " << error.what();
      if(modulesConfig.Has<DelphesParameters>(moduleName))
        message << "\nParameters:\n"
                << modulesConfig.Get<DelphesParameters>(moduleName);
      throw std::runtime_error(message.str());
    }
  }
}

//------------------------------------------------------------------------------

void PyDelphes::SetReaderConfig(const pybind11::dict &readerArgs)
{
  fReaderConfig = readerArgs;
  const DelphesParameters readerConfig = fReaderConfig.cast<DelphesParameters>();
  fEventReader = DelphesReaderFactory::Get().Build(readerConfig.Get<std::string>("ReaderType"), readerConfig);
  if(readerArgs.contains("inputFiles"))
  {
    if(const std::vector<std::string> inputFiles = readerArgs["inputFiles"].cast<std::vector<std::string> >(); !inputFiles.empty())
      fEventReader->LoadInputFile(inputFiles.at(0));
    //TODO: manage multi-files inputs
  }
  fEventReader->SetFactory(GetFactory());
}

//------------------------------------------------------------------------------

void PyDelphes::LoadTCL(std::string_view tclFilePath)
{
  DelphesTCLConfReader confReader;
  confReader.ReadFile(tclFilePath);
  fModulesConfig = py::cast(confReader.Parameters());
}

//------------------------------------------------------------------------------

const PyDelphesEvent &PyDelphes::Next()
{
  if(!fEventReader)
    throw std::runtime_error("Event reader object is not yet initialised.");
  if(!fIsInitialised)
  {
    InitTask();
    fLastModulesConfig = fModulesConfig;
    fIsInitialised = true;
  }
  Clear();
  fEventReader->Clear();
  fEventReader->ReadEvent();
  ProcessTask();
  return *fDelphesEvent;
}

//------------------------------------------------------------------------------

void PyDelphes::SetModulesConfig(const py::dict &modulesObj)
{
  py::list execPath;
  for(const std::pair<py::handle, py::handle> &moduleObj : modulesObj)
  {
    const std::string moduleName = moduleObj.first.cast<std::string>();
    if(py::isinstance<py::dict>(moduleObj.second))
    {
      fModulesConfig[py::str(moduleName)] = moduleObj.second.cast<py::dict>();
      execPath.append(py::str(moduleName));
    }
  }
  fModulesConfig[py::str("ExecutionPath")] = execPath;
  fIsInitialised = false;
}

//------------------------------------------------------------------------------
