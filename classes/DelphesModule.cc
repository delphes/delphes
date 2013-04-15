
/** \class DelphesModule
 *
 *  Base class for all Delphes modules
 *
 *  $Date: 2008-06-04 13:57:25 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include "classes/DelphesFactory.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TROOT.h"
#include "TClass.h"
#include "TFolder.h"
#include "TObjArray.h"

#include <iostream>
#include <stdexcept>
#include <sstream>

using namespace std;

DelphesModule::DelphesModule() :
  fTreeWriter(0), fFactory(0), fPlots(0),
  fPlotFolder(0), fExportFolder(0)
{
}

//------------------------------------------------------------------------------

DelphesModule::~DelphesModule()
{
}

//------------------------------------------------------------------------------

void DelphesModule::Init()
{
}

//------------------------------------------------------------------------------

void DelphesModule::Process()
{
}

//------------------------------------------------------------------------------

void DelphesModule::Finish()
{
}

//------------------------------------------------------------------------------

TObjArray *DelphesModule::ImportArray(const char *name)
{
  stringstream message;
  TObjArray *object;

  object = static_cast<TObjArray *>(GetObject(Form("Export/%s", name), TObjArray::Class()));
  if(!object)
  {
    message << "can't access input list '" << name;
    message << "' in module '" << GetName() << "'";
    throw runtime_error(message.str());
  }

  return object;
}

//------------------------------------------------------------------------------

TObjArray *DelphesModule::ExportArray(const char *name)
{
  TObjArray *array;
  if(!fExportFolder)
  {
    fExportFolder = NewFolder("Export");
  }

  array = GetFactory()->NewPermanentArray();

  array->SetName(name);
  fExportFolder->Add(array);

  return array;
}

//------------------------------------------------------------------------------

ExRootTreeBranch *DelphesModule::NewBranch(const char *name, TClass *cl)
{
  stringstream message;
  if(!fTreeWriter)
  {
    fTreeWriter = static_cast<ExRootTreeWriter *>(GetObject("TreeWriter", ExRootTreeWriter::Class()));
    if(!fTreeWriter)
    {
      message << "can't access access tree writer";
      throw runtime_error(message.str());
    }
  }
  return fTreeWriter->NewBranch(name, cl);
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

DelphesFactory *DelphesModule::GetFactory()
{
  stringstream message;
  if(!fFactory)
  {
    fFactory = static_cast<DelphesFactory *>(GetObject("ObjectFactory", DelphesFactory::Class()));
    if(!fFactory)
    {
      message << "can't access access object factory";
      throw runtime_error(message.str());
    }
  }
  return fFactory;
}


