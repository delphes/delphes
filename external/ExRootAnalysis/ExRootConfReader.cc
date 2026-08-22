
/** \class ExRootConfReader
 *
 *  Class handling configuration data
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootConfReader.h"

#include "tcl/jim.h"

#include "TSystem.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace std;

static Jim_CmdProc SourceCmdProc;

//------------------------------------------------------------------------------

ExRootConfReader::ExRootConfReader() :
  fTopDir(0), fTclInterp(0)
{
  fTclInterp = Jim_CreateInterp();
  Jim_RegisterCoreCommands(fTclInterp);
  Jim_InitStaticExtensions(fTclInterp);

  Jim_CreateCommand(fTclInterp, "source", SourceCmdProc, this, NULL);
}

//------------------------------------------------------------------------------

ExRootConfReader::~ExRootConfReader()
{
  Jim_FreeInterp(fTclInterp);
}

//------------------------------------------------------------------------------

void ExRootConfReader::ReadData(const char *dirName, char *data, int length)
{
  stringstream message;
  Jim_Obj *object;
  int rc;

  fTopDir = dirName;

  object = Jim_NewStringObj(fTclInterp, data, length);
  object->bytes = data;
  object->length = length;
  object->typePtr = NULL;

  Jim_IncrRefCount(object);
  rc = Jim_EvalObj(fTclInterp, object);
  object->bytes = NULL;
  object->length = 0;
  Jim_DecrRefCount(fTclInterp, object);

  if(rc != JIM_OK)
  {
    message << "can't read configuration data" << endl;
    if(rc == JIM_BREAK)
    {
      message << "invoked \"break\" outside of a loop" << endl;
    }
    else if(rc == JIM_CONTINUE)
    {
      message << "invoked \"continue\" outside of a loop" << endl;
    }
    else
    {
      Jim_MakeErrorMessage(fTclInterp);
      message << Jim_GetString(Jim_GetResult(fTclInterp), NULL);
    }
    throw runtime_error(message.str());
  }
}

//------------------------------------------------------------------------------

void ExRootConfReader::ReadFile(const char *fileName, bool isTop)
{
  stringstream message;
  Jim_Obj *object;
  int length, rc;
  char *buffer;

  ifstream inputFileStream(fileName, ios::in | ios::ate);
  if(!inputFileStream.is_open())
  {
    message << "can't open configuration file " << fileName;
    throw runtime_error(message.str());
  }

  if(isTop) fTopDir = gSystem->DirName(fileName);

  length = inputFileStream.tellg();
  inputFileStream.seekg(0, ios::beg);
  inputFileStream.clear();
  buffer = new char[length + 1];
  buffer[length] = 0;
  inputFileStream.read(buffer, length);

  object = Jim_NewObj(fTclInterp);
  object->bytes = buffer;
  object->length = length;
  object->typePtr = NULL;

  Jim_IncrRefCount(object);
  rc = Jim_EvalObj(fTclInterp, object);
  object->bytes = NULL;
  object->length = 0;
  Jim_DecrRefCount(fTclInterp, object);

  delete[] buffer;

  if(rc != JIM_OK)
  {
    message << "can't read configuration file " << fileName << endl;
    if(rc == JIM_BREAK)
    {
      message << "invoked \"break\" outside of a loop" << endl;
    }
    else if(rc == JIM_CONTINUE)
    {
      message << "invoked \"continue\" outside of a loop" << endl;
    }
    else
    {
      Jim_MakeErrorMessage(fTclInterp);
      message << Jim_GetString(Jim_GetResult(fTclInterp), NULL);
    }
    throw runtime_error(message.str());
  }
}

//------------------------------------------------------------------------------

ExRootConfParam ExRootConfReader::GetParam(const char *name)
{
  Jim_Obj *object;
  Jim_Obj *variableName = Jim_NewStringObj(fTclInterp, const_cast<char *>(name), -1);
  object = Jim_GetGlobalVariable(fTclInterp, variableName, 0);
  return ExRootConfParam(name, object, fTclInterp);
}

//------------------------------------------------------------------------------

int ExRootConfReader::GetInt(const char *name, int defaultValue, int index)
{
  ExRootConfParam object = GetParam(name);
  if(index >= 0)
  {
    object = object[index];
  }

  return object.GetInt(defaultValue);
}

//------------------------------------------------------------------------------

long ExRootConfReader::GetLong(const char *name, long defaultValue, int index)
{
  ExRootConfParam object = GetParam(name);
  if(index >= 0)
  {
    object = object[index];
  }

  return object.GetLong(defaultValue);
}

//------------------------------------------------------------------------------

double ExRootConfReader::GetDouble(const char *name, double defaultValue, int index)
{
  ExRootConfParam object = GetParam(name);
  if(index >= 0)
  {
    object = object[index];
  }

  return object.GetDouble(defaultValue);
}

//------------------------------------------------------------------------------

bool ExRootConfReader::GetBool(const char *name, bool defaultValue, int index)
{
  ExRootConfParam object = GetParam(name);
  if(index >= 0)
  {
    object = object[index];
  }

  return object.GetBool(defaultValue);
}

//------------------------------------------------------------------------------

const char *ExRootConfReader::GetString(const char *name, const char *defaultValue, int index)
{
  ExRootConfParam object = GetParam(name);
  if(index >= 0)
  {
    object = object[index];
  }

  return object.GetString(defaultValue);
}

//------------------------------------------------------------------------------

int SourceCmdProc(Jim_Interp *interp, int objc, Jim_Obj *const objv[])
{
  ExRootConfReader *reader = static_cast<ExRootConfReader *>(Jim_CmdPrivData(interp));
  stringstream fileName;

  if(objc != 2)
  {
    Jim_WrongNumArgs(interp, 1, objv, "fileName");
    return JIM_ERR;
  }

  fileName << reader->GetTopDir() << "/" << Jim_GetString(objv[1], NULL);

  try
  {
    reader->ReadFile(fileName.str().c_str(), false);
  }
  catch(runtime_error &e)
  {
    Jim_SetResultString(interp, e.what(), -1);
    return JIM_ERR;
  }

  return JIM_OK;
}

//------------------------------------------------------------------------------

ExRootConfParam::ExRootConfParam(const char *name, Jim_Obj *object, Jim_Interp *interp) :
  TNamed(name, ""), fObject(object), fTclInterp(interp)
{
}

//------------------------------------------------------------------------------

int ExRootConfParam::GetInt(int defaultValue)
{
  stringstream message;
  long result;
  int rc;

  if(!fObject) return defaultValue;

  rc = Jim_GetLong(fTclInterp, fObject, &result);

  if(rc == JIM_OK && static_cast<int>(result) == result) return result;

  message << "parameter '" << GetName() << "' is not an integer." << endl;
  message << GetName() << " = " << Jim_GetString(fObject, NULL);
  throw runtime_error(message.str());
}

//------------------------------------------------------------------------------

long ExRootConfParam::GetLong(long defaultValue)
{
  stringstream message;
  long result;
  int rc;

  if(!fObject) return defaultValue;

  rc = Jim_GetLong(fTclInterp, fObject, &result);

  if(rc == JIM_OK) return result;

  message << "parameter '" << GetName() << "' is not an long integer." << endl;
  message << GetName() << " = " << Jim_GetString(fObject, NULL);
  throw runtime_error(message.str());
}

//------------------------------------------------------------------------------

double ExRootConfParam::GetDouble(double defaultValue)
{
  stringstream message;
  double result;
  int rc;

  if(!fObject) return defaultValue;

  rc = Jim_GetDouble(fTclInterp, fObject, &result);

  if(rc == JIM_OK) return result;

  message << "parameter '" << GetName() << "' is not a number." << endl;
  message << GetName() << " = " << Jim_GetString(fObject, NULL);
  throw runtime_error(message.str());
}

//------------------------------------------------------------------------------

bool ExRootConfParam::GetBool(bool defaultValue)
{
  stringstream message;
  int result;
  int rc;

  if(!fObject) return defaultValue;

  rc = Jim_GetBoolean(fTclInterp, fObject, &result);

  if(rc == JIM_OK) return result;

  message << "parameter '" << GetName() << "' is not a boolean." << endl;
  message << GetName() << " = " << Jim_GetString(fObject, NULL);
  throw runtime_error(message.str());
}

//------------------------------------------------------------------------------

const char *ExRootConfParam::GetString(const char *defaultValue)
{
  if(!fObject) return defaultValue;

  return Jim_GetString(fObject, NULL);
}

//------------------------------------------------------------------------------

int ExRootConfParam::GetSize()
{
  if(!fObject) return 0;

  return Jim_ListLength(fTclInterp, fObject);
}

//------------------------------------------------------------------------------

ExRootConfParam ExRootConfParam::operator[](int index)
{
  stringstream message;
  Jim_Obj *object;
  int rc;

  if(!fObject) return ExRootConfParam(GetName(), 0, fTclInterp);

  rc = Jim_ListIndex(fTclInterp, fObject, index, &object, 0);

  if(rc == JIM_OK) return ExRootConfParam(GetName(), object, fTclInterp);

  message << "list index for parameter '" << GetName() << "' is out of range." << endl;
  message << GetName() << " = " << Jim_GetString(fObject, NULL);
  throw runtime_error(message.str());
}
