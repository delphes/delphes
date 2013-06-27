
/** \class Delphes
 *
 *  Main Delphes module.
 *  Controls execution of all other modules.
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Delphes.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootConfReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "TROOT.h"
#include "TMath.h"
#include "TFolder.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm> 
#include <stdexcept>
#include <iostream>
#include <sstream>

#include <string.h>
#include <stdio.h>

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
  if(fFactory) delete fFactory;
  TFolder *folder = GetFolder();
  if(folder)
  {
    folder->Clear();
    delete folder;
  }
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
