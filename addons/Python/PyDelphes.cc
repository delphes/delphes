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

#include "classes/DelphesFactory.h"
#include "classes/DelphesReader.h"
#include "modules/Delphes.h"

#include "PyDelphes.h"
#include "PyDelphesConfReader.h"

#include <ExRootAnalysis/ExRootProgressBar.h>

namespace py = pybind11;

PyDelphes::~PyDelphes() { FinishTask(); }

//------------------------------------------------------------------------------

void PyDelphes::Init()
{
  for(const std::string &moduleName : fConfig.Get<std::vector<std::string> >("ExecutionPath"))
  {
    if(!fConfig.Has<DelphesParameters>(moduleName))
    {
      std::ostringstream message;
      message << "module '" << moduleName;
      message << "' is specified in ExecutionPath but not configured.";
      throw std::runtime_error(message.str());
    }
    DelphesFactory *factory = GetFactory();
    if(!factory)
      throw std::runtime_error("Delphes factory was not initialised for this PyDelphes module.");
    fDelphesEvent = std::make_unique<PyDelphesEvent>(*factory);
    try
    {
      const DelphesParameters moduleParams = fConfig.Get<DelphesParameters>(moduleName);
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
      if(fConfig.Has<DelphesParameters>(moduleName))
        message << "\nParameters:\n"
                << fConfig.Get<DelphesParameters>(moduleName);
      throw std::runtime_error(message.str());
    }
  }
}

//------------------------------------------------------------------------------

void PyDelphes::SetReader(DelphesReader *readerObj)
{
  Delphes::SetReader(readerObj);
  readerObj->SetFactory(GetFactory());
  readerObj->SetMaxEvents(fConfig.Get<int>("MaxEvents", 0));
  readerObj->SetSkipEvents(fConfig.Get<int>("SkipEvents", 0));
  readerObj->Clear();
}

//------------------------------------------------------------------------------

const PyDelphesEvent &PyDelphes::Next()
{
  DelphesReader *eventReader = GetReader();
  if(!eventReader)
    throw std::runtime_error("Event reader object is not yet initialised.");
  if(!fIsInitialised)
  {
    InitTask();
    fIsInitialised = true;
  }
  Clear();
  std::cout << __PRETTY_FUNCTION__ << ":" << eventReader << std::endl;
  eventReader->ReadEvent();
  ProcessTask();
  Clear();
  eventReader->Clear();
  return *fDelphesEvent;
}

//------------------------------------------------------------------------------

py::dict PyDelphes::GetModules() const { return py::dict{}; /*fConfig*/ } //TODO: make UI a bit better

//------------------------------------------------------------------------------

void PyDelphes::SetModules(const py::dict &modulesObj)
{
  for(const std::pair<py::handle, py::handle> &moduleObj : modulesObj)
  {
    const std::string moduleName = moduleObj.first.cast<std::string>();
    if(py::isinstance<py::dict>(moduleObj.second))
      fConfig.Set(moduleName, PyDelphesConfReader{moduleObj.second.cast<py::dict>()}.Parameters());
  }
}

//------------------------------------------------------------------------------
