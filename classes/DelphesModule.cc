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

/** \class DelphesModule
 *
 *  Base class for all Delphes modules
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesModule.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "TClass.h"
#include "TFolder.h"

#include <sstream>
#include <stdexcept>

using namespace std;

CandidatesCollection DelphesModule::ImportArray(const char *name)
{
  DelphesFactory *factory = GetFactory();
  if(!factory)
    throw std::runtime_error("Failed to retrieve the Delphes objects factory for module '" + std::string{GetName()} + "'.");
  if(!factory->Has(name))
  {
    std::ostringstream message;
    message << "can't access input list '" << name;
    message << "' in module '" << GetName() << "'";
    throw std::runtime_error(message.str());
  }
  return factory->Attach<std::vector<Candidate *> >(name);
}

//------------------------------------------------------------------------------

CandidatesCollection DelphesModule::ExportArray(const char *name)
{
  DelphesFactory *factory = GetFactory();
  if(!factory)
    throw std::runtime_error("Failed to retrieve the Delphes objects factory for module '" + std::string{GetName()} + "'.");
  std::ostringstream collectionLabel;
  collectionLabel << GetName() << "/" << name;
  if(factory->Has(name))
  {
    std::ostringstream message;
    message << "Collection with name '" << name << "' and label '" << collectionLabel.str() << "'  was already booked in the memory slot.";
    throw std::runtime_error(message.str());
  }
  return factory->Book<std::vector<Candidate *> >(collectionLabel.str());
}

//------------------------------------------------------------------------------

ExRootTreeBranch *DelphesModule::NewBranch(const char *name, TClass *cl)
{
  return GetTreeWriter()->NewBranch(name, cl);
}

//------------------------------------------------------------------------------

void DelphesModule::AddInfo(const char *name, Double_t value)
{
  GetTreeWriter()->AddInfo(name, value);
}

//------------------------------------------------------------------------------

ExRootResult *DelphesModule::GetPlots()
{
  if(!fPlots)
  {
    fPlots = new ExRootResult();
    fPlots->SetFolder(GetFolder());
  }
  return fPlots;
}

//------------------------------------------------------------------------------

DelphesFactory *DelphesModule::GetFactory() const
{
  if(!fFactory)
  {
    std::ostringstream message;
    message << "can't access access object factory for module '" << GetName() << "'.";
    throw std::runtime_error(message.str());
  }
  return fFactory;
}

//------------------------------------------------------------------------------

ExRootTreeWriter *DelphesModule::GetTreeWriter() const
{
  if(!fTreeWriter)
  {
    std::ostringstream message;
    message << "can't access access tree writer for module '" << GetName() << "'.";
    throw std::runtime_error(message.str());
  }
  return fTreeWriter;
}

//------------------------------------------------------------------------------
