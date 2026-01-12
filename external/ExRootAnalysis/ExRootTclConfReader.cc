
/** \class ExRootConfReader
 *
 *  Class handling output ROOT tree
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootTclConfReader.h"

#include "tcl/tcl.h"

#include "TSystem.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace std;

static Tcl_ObjCmdProc ModuleObjCmdProc;
static Tcl_ObjCmdProc SourceObjCmdProc;

//------------------------------------------------------------------------------

ExRootTclConfReader::ExRootTclConfReader() :
  fTclInterp(Tcl_CreateInterp())
{
  Tcl_CreateObjCommand(fTclInterp, const_cast<char *>("module"), ModuleObjCmdProc, this, 0);
  Tcl_CreateObjCommand(fTclInterp, const_cast<char *>("source"), SourceObjCmdProc, this, 0);
}

//------------------------------------------------------------------------------

ExRootTclConfReader::~ExRootTclConfReader()
{
  Tcl_DeleteInterp(fTclInterp);
}

//------------------------------------------------------------------------------

void ExRootTclConfReader::ReadFile(const char *fileName, bool isTop)
{
  stringstream message;
  int length;
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

  Tcl_Obj *cmdObjPtr = Tcl_NewObj();
  cmdObjPtr->bytes = buffer;
  cmdObjPtr->length = length;

  Tcl_IncrRefCount(cmdObjPtr);

  if(Tcl_EvalObj(fTclInterp, cmdObjPtr) != TCL_OK)
  {
    message << "can't read configuration file " << fileName << endl;
    message << Tcl_GetStringResult(fTclInterp);
    throw runtime_error(message.str());
  }

  cmdObjPtr->bytes = 0;
  cmdObjPtr->length = 0;

  Tcl_DecrRefCount(cmdObjPtr);

  delete[] buffer;
}

//------------------------------------------------------------------------------

std::unique_ptr<ExRootConfParam> ExRootTclConfReader::GetParam(const char *name)
{
  Tcl_Obj *object;
  Tcl_Obj *variableName = Tcl_NewStringObj(const_cast<char *>(name), -1);
  object = Tcl_ObjGetVar2(fTclInterp, variableName, 0, TCL_GLOBAL_ONLY);
  return std::make_unique<ExRootTclConfParam>(name, object, fTclInterp);
}

//------------------------------------------------------------------------------

std::unique_ptr<ExRootConfParam> ExRootTclConfReader::GetGlobalParam(const char *name) { return GetParam(TString("::") + name); }

//------------------------------------------------------------------------------

int ModuleObjCmdProc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[])
{
  ExRootTclConfReader *reader = static_cast<ExRootTclConfReader *>(clientData);

  if(objc < 3)
  {
    Tcl_WrongNumArgs(interp, 1, objv, "className moduleName ?arg...?");
    return TCL_ERROR;
  }

  // add module to a list of modules to be created

  reader->AddModule(Tcl_GetStringFromObj(objv[1], 0), Tcl_GetStringFromObj(objv[2], 0));

  if(objc > 3)
  {
    Tcl_Obj *object = Tcl_NewListObj(0, 0);
    Tcl_ListObjAppendElement(interp, object, Tcl_NewStringObj("namespace", -1));
    Tcl_ListObjAppendElement(interp, object, Tcl_NewStringObj("eval", -1));
    Tcl_ListObjAppendList(interp, object, Tcl_NewListObj(objc - 2, objv + 2));

    return Tcl_GlobalEvalObj(interp, object);
  }

  return TCL_OK;
}

//------------------------------------------------------------------------------

int SourceObjCmdProc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[])
{
  ExRootTclConfReader *reader = static_cast<ExRootTclConfReader *>(clientData);
  stringstream fileName;

  if(objc != 2)
  {
    Tcl_WrongNumArgs(interp, 1, objv, "fileName");
    return TCL_ERROR;
  }

  fileName << reader->GetTopDir() << "/" << Tcl_GetStringFromObj(objv[1], 0);
  reader->ReadFile(fileName.str().c_str(), false);

  return TCL_OK;
}

//------------------------------------------------------------------------------

ExRootTclConfParam::ExRootTclConfParam(const char *name, Tcl_Obj *object, Tcl_Interp *interp) :
  ExRootConfParam(name), fObject(object), fTclInterp(interp)
{
}

//------------------------------------------------------------------------------

int ExRootTclConfParam::GetInt(int defaultValue)
{
  stringstream message;
  int result = defaultValue;
  if(fObject && TCL_OK != Tcl_GetIntFromObj(fTclInterp, fObject, &result))
  {
    message << "parameter '" << fName << "' is not an integer." << endl;
    message << fName << " = " << Tcl_GetStringFromObj(fObject, 0);
    throw runtime_error(message.str());
  }
  return result;
}

//------------------------------------------------------------------------------

long ExRootTclConfParam::GetLong(long defaultValue)
{
  stringstream message;
  long result = defaultValue;
  if(fObject && TCL_OK != Tcl_GetLongFromObj(fTclInterp, fObject, &result))
  {
    message << "parameter '" << fName << "' is not an long integer." << endl;
    message << fName << " = " << Tcl_GetStringFromObj(fObject, 0);
    throw runtime_error(message.str());
  }
  return result;
}

//------------------------------------------------------------------------------

double ExRootTclConfParam::GetDouble(double defaultValue)
{
  stringstream message;
  double result = defaultValue;
  if(fObject && TCL_OK != Tcl_GetDoubleFromObj(fTclInterp, fObject, &result))
  {
    message << "parameter '" << fName << "' is not a number." << endl;
    message << fName << " = " << Tcl_GetStringFromObj(fObject, 0);
    throw runtime_error(message.str());
  }
  return result;
}

//------------------------------------------------------------------------------

bool ExRootTclConfParam::GetBool(bool defaultValue)
{
  stringstream message;
  int result = defaultValue;
  if(fObject && TCL_OK != Tcl_GetBooleanFromObj(fTclInterp, fObject, &result))
  {
    message << "parameter '" << fName << "' is not a boolean." << endl;
    message << fName << " = " << Tcl_GetStringFromObj(fObject, 0);
    throw runtime_error(message.str());
  }
  return result;
}

//------------------------------------------------------------------------------

const char *ExRootTclConfParam::GetString(const char *defaultValue)
{
  const char *result = defaultValue;
  if(fObject) result = Tcl_GetStringFromObj(fObject, 0);
  return result;
}

//------------------------------------------------------------------------------

int ExRootTclConfParam::GetSize()
{
  stringstream message;
  int length = 0;
  if(fObject && TCL_OK != Tcl_ListObjLength(fTclInterp, fObject, &length))
  {
    message << "parameter '" << fName << "' is not a list." << endl;
    message << fName << " = " << Tcl_GetStringFromObj(fObject, 0);
    throw runtime_error(message.str());
  }
  return length;
}

//------------------------------------------------------------------------------

std::unique_ptr<ExRootConfParam> ExRootTclConfParam::operator[](int index)
{
  Tcl_Obj *object = 0;
  if(fObject && TCL_OK != Tcl_ListObjIndex(fTclInterp, fObject, index, &object))
  {
    stringstream message;
    message << "parameter '" << fName << "' is not a list." << endl;
    message << fName << " = " << Tcl_GetStringFromObj(fObject, 0);
    throw runtime_error(message.str());
  }
  return std::make_unique<ExRootTclConfParam>(fName, object, fTclInterp);
}

std::unique_ptr<ExRootConfParam> ExRootTclConfParam::GetParam(const char *paramName) const
{
  //Tcl_Obj *object = 0;
  //if(fObject && TCL_OK != Tcl_DictObjGet(fTclInterp, fObject, paramName, &object))
  //TODO: bump version of TCL to enable dictionary objects parsing
  {
    stringstream message;
    message << "parameter '" << fName << "' is not a dictionary." << endl;
    message << fName << " = " << Tcl_GetStringFromObj(fObject, 0);
    throw runtime_error(message.str());
  }
  //return std::make_unique<ExRootTclConfParam>(fName, object, fTclInterp);
}
