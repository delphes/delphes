
/** \class ExRootConfReader
 *
 *  Class handling input steering card
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootPythonConfReader.h"
#include "ExRootAnalysis/ExRootTclConfReader.h"

#include "TSystem.h"

#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

//------------------------------------------------------------------------------

ExRootConfReader::ExRootConfReader() :
  fTopDir(0)
{
}

//------------------------------------------------------------------------------

std::unique_ptr<ExRootConfReader> ExRootConfReader::ReadConf(const char *fileName)
{
  if(const size_t str_len = strlen(fileName); // first check if we have a .tcl extension
    str_len >= strlen(".tcl") && (0 == strcmp(fileName + (str_len - strlen(".tcl")), ".tcl")))
  {
    auto reader = std::make_unique<ExRootTclConfReader>();
    reader->ReadFile(fileName);
    return reader;
  }
  else if(str_len >= strlen(".py") && (0 == strcmp(fileName + (str_len - strlen(".py")), ".py")))
  {
    auto reader = std::make_unique<ExRootPythonConfReader>();
    reader->ReadFile(fileName);
    return reader;
  }
  ostringstream message;
  message << "Unsupported extension for input configuration card '" << fileName << "'.";
  throw runtime_error(message.str());
}

//------------------------------------------------------------------------------

int ExRootConfReader::GetInt(const char *name, int defaultValue, int index)
{
  auto object = GetParam(name);
  if(index >= 0)
    object = (*object)[index];
  return object->GetInt(defaultValue);
}

//------------------------------------------------------------------------------

long ExRootConfReader::GetLong(const char *name, long defaultValue, int index)
{
  auto object = GetParam(name);
  if(index >= 0)
    object = (*object)[index];
  return object->GetLong(defaultValue);
}

//------------------------------------------------------------------------------

double ExRootConfReader::GetDouble(const char *name, double defaultValue, int index)
{
  auto object = GetParam(name);
  if(index >= 0)
    object = (*object)[index];
  return object->GetDouble(defaultValue);
}

//------------------------------------------------------------------------------

bool ExRootConfReader::GetBool(const char *name, bool defaultValue, int index)
{
  auto object = GetParam(name);
  if(index >= 0)
    object = (*object)[index];
  return object->GetBool(defaultValue);
}

//------------------------------------------------------------------------------

const char *ExRootConfReader::GetString(const char *name, const char *defaultValue, int index)
{
  auto object = GetParam(name);
  if(index >= 0)
    object = (*object)[index];
  return object->GetString(defaultValue);
}

//------------------------------------------------------------------------------

std::unique_ptr<ExRootConfParam> ExRootConfReader::GetGlobalParam(const char *name) { return GetParam(name); }

//------------------------------------------------------------------------------

void ExRootConfReader::AddModule(const char *className, const char *moduleName)
{
  if(ExRootTaskMap::iterator itModules = fModules.find(moduleName); itModules != fModules.end())
  {
    cout << "** WARNING: module '" << moduleName << "' is already configured.";
    cout << " Only first entry will be used." << endl;
  }
  else
  {
    fModules.insert(make_pair(moduleName, className));
    cout << left;
    cout << setw(30) << "** INFO: adding module";
    cout << setw(25) << className;
    cout << setw(25) << moduleName << endl;
  }
}

//------------------------------------------------------------------------------

ExRootConfParam::ExRootConfParam(const char *name) :
  fName(name)
{
}
