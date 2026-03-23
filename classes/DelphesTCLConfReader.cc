/** \class DelphesTCLConfReader
 *
 *  Class handling TCL input card parsing
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesTCLConfReader.h"

#include <tcl/tcl.h>

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

DelphesTCLConfReader::DelphesTCLConfReader() :
  fTclInterp(Tcl_CreateInterp())
{
  Tcl_CreateObjCommand(fTclInterp.get(), const_cast<char *>("module"), ModuleObjCmdProc, this, 0);
  Tcl_CreateObjCommand(fTclInterp.get(), const_cast<char *>("source"), SourceObjCmdProc, this, 0);
}

//------------------------------------------------------------------------------

std::string DelphesTCLConfReader::Run(std::string_view command) const
{
  if(Tcl_Eval(fTclInterp.get(), const_cast<char *>(command.data())) != TCL_OK) return "nil";
  return Tcl_GetStringFromObj(Tcl_GetObjResult(fTclInterp.get()), 0);
}

//------------------------------------------------------------------------------

void DelphesTCLConfReader::ReadFile(std::string_view fileName)
{
  std::ifstream inputFileStream(std::string{fileName}, ios::in | ios::ate);
  if(!inputFileStream.is_open())
  {
    std::ostringstream message;
    message << "can't open configuration file " << fileName;
    throw std::runtime_error(message.str());
  }

  const int length = inputFileStream.tellg();
  inputFileStream.seekg(0, ios::beg);
  inputFileStream.clear();
  std::vector<char> buffer(length + 1, 0);
  inputFileStream.read(buffer.data(), length);

  Tcl_Obj *cmdObjPtr = Tcl_NewObj();
  cmdObjPtr->bytes = buffer.data();
  cmdObjPtr->length = length;

  Tcl_IncrRefCount(cmdObjPtr);

  if(Tcl_EvalObj(fTclInterp.get(), cmdObjPtr) != TCL_OK)
  {
    std::ostringstream message;
    message << "can't read configuration file " << fileName << std::endl;
    message << Tcl_GetStringResult(fTclInterp.get());
    throw std::runtime_error(message.str());
  }

  cmdObjPtr->bytes = 0;
  cmdObjPtr->length = 0;

  Tcl_DecrRefCount(cmdObjPtr);

  ParseValue(nullptr, "::ExecutionPath", "ExecutionPath", fParams);

  for(const std::string &moduleName : fParams.Get<std::vector<std::string> >("ExecutionPath"))
  {
    DelphesParameters moduleParams;
    ParseParameters(moduleName, moduleParams);
    moduleParams.Set("ModuleType", fModuleTypes.at(moduleName));
    fParams.Set(moduleName, moduleParams);
  }

  //std::cout << ">>> " << Run("namespace children ::ScalarHT") << std::endl;
  //std::cout << ">>> " << Run("info vars ::ScalarHT::*") << std::endl;

  //std::cout << fParams << std::endl;
}

//------------------------------------------------------------------------------

void DelphesTCLConfReader::ParseValue(const Tcl_Obj *tclObject,
  std::string_view tclName, std::string_view keyName,
  DelphesParameters &delphesParams) const
{
  Tcl_Obj *object = nullptr;
  if(tclObject)
    object = const_cast<Tcl_Obj *>(tclObject);
  else if(int length; tclObject && Tcl_ListObjLength(fTclInterp.get(), const_cast<Tcl_Obj *>(tclObject), &length) == TCL_OK) // we have a list/collection
  {
    throw std::runtime_error("Dictionary-like browsing is not supported by this version of TCL.");
  }
  else // retrieve from the global environment
  {
    Tcl_Obj *variableName = Tcl_NewStringObj(const_cast<char *>(tclName.data()), -1); // build the key name object
    object = Tcl_ObjGetVar2(fTclInterp.get(), variableName, 0, TCL_GLOBAL_ONLY);
  }
  if(!object)
  {
    std::ostringstream message;
    message << "Failed to retrieve object with key '" << tclName << "' from object.";
    throw std::runtime_error(message.str());
  }

  const std::string delphesKeyName = keyName.empty() ? std::string{tclName} : std::string{keyName};

  if(int length; Tcl_ListObjLength(fTclInterp.get(), object, &length) == TCL_OK)
  { // we have a list/collection
    // with this version of TCL, only flat lists are supported, no need to parse dictionaries
    std::vector<long> longCollection;
    std::vector<int> intCollection;
    std::vector<double> doubleCollection;
    std::vector<bool> boolCollection;
    std::vector<std::string> stringCollection;
    for(int i = 0; i < length; ++i)
    {
      if(Tcl_Obj *subObject = nullptr; Tcl_ListObjIndex(fTclInterp.get(), object, i, &subObject) == TCL_OK)
      {
        if(long result; Tcl_GetLongFromObj(fTclInterp.get(), subObject, &result) == TCL_OK)
          longCollection.emplace_back(result);
        else if(int result; Tcl_GetIntFromObj(fTclInterp.get(), subObject, &result) == TCL_OK)
          intCollection.emplace_back(result);
        else if(double result; Tcl_GetDoubleFromObj(fTclInterp.get(), subObject, &result) == TCL_OK)
          doubleCollection.emplace_back(result);
        else if(int result; Tcl_GetBooleanFromObj(fTclInterp.get(), object, &result) == TCL_OK)
          boolCollection.emplace_back(result);
        else if(std::string result = Tcl_GetStringFromObj(subObject, 0); !result.empty())
          stringCollection.emplace_back(result);
      }
    }
    if(!longCollection.empty()) delphesParams.Set(delphesKeyName, longCollection);
    if(!intCollection.empty()) delphesParams.Set(delphesKeyName, intCollection);
    if(!doubleCollection.empty()) delphesParams.Set(delphesKeyName, doubleCollection);
    if(!boolCollection.empty()) delphesParams.Set(delphesKeyName, boolCollection);
    if(!stringCollection.empty())
    {
      if(delphesKeyName.find("Formula") != std::string::npos)
      { //FIXME: hacky way to ensure formulas are parsed as full strings
        if(std::string result = Tcl_GetStringFromObj(object, 0); !result.empty())
        { //FIXME check if reasonable
          delphesParams.Set(delphesKeyName, result);
          return;
        }
      }
      delphesParams.Set(delphesKeyName, stringCollection);
    }
  }
  else if(long result; Tcl_GetLongFromObj(fTclInterp.get(), object, &result) == TCL_OK)
    delphesParams.Set(delphesKeyName, result);
  else if(int result; Tcl_GetIntFromObj(fTclInterp.get(), object, &result) == TCL_OK)
    delphesParams.Set(delphesKeyName, result);
  else if(double result; Tcl_GetDoubleFromObj(fTclInterp.get(), object, &result) == TCL_OK)
    delphesParams.Set(delphesKeyName, result);
  else if(int result; Tcl_GetBooleanFromObj(fTclInterp.get(), object, &result) == TCL_OK)
    delphesParams.Set(delphesKeyName, static_cast<bool>(result));
  else if(std::string result = Tcl_GetStringFromObj(object, 0); !result.empty()) //FIXME check if reasonable
    delphesParams.Set(delphesKeyName, result);
}

//------------------------------------------------------------------------------

std::vector<std::string> DelphesTCLConfReader::ChildrenList(std::string_view namespaceName) const
{
  std::ostringstream command;
  //command << "namespace children " << namespaceName;
  command << "info vars " << namespaceName << "::*";
  if(Tcl_Eval(fTclInterp.get(), const_cast<char *>(command.str().data())) != TCL_OK) return {};
  Tcl_Obj *namespaceChildren = Tcl_GetObjResult(fTclInterp.get());
  int length;
  if(Tcl_ListObjLength(fTclInterp.get(), namespaceChildren, &length) != TCL_OK)
  {
    std::ostringstream message;
    message << "Returned object for '" << namespaceName << "' namespace children evaluation is not a collection.";
    std::runtime_error(message.str());
  }
  std::vector<std::string> childrenList(length);
  for(int i = 0; i < length; ++i)
  {
    if(Tcl_Obj *subObject = nullptr; Tcl_ListObjIndex(fTclInterp.get(), namespaceChildren, i, &subObject) == TCL_OK && subObject)
      childrenList[i] = Tcl_GetStringFromObj(subObject, 0);
  }
  return childrenList;
}

//------------------------------------------------------------------------------

std::string DelphesTCLConfReader::TrimmedName(std::string_view keyName) const
{
  std::ostringstream command;
  command << "namespace tail " << keyName;
  if(Tcl_Eval(fTclInterp.get(), const_cast<char *>(command.str().data())) == TCL_OK)
  {
    if(Tcl_Obj *trimmedSubObject = Tcl_GetObjResult(fTclInterp.get()); trimmedSubObject)
      return Tcl_GetStringFromObj(trimmedSubObject, 0);
  }
  return std::string{keyName};
}

//------------------------------------------------------------------------------

void DelphesTCLConfReader::ParseParameters(std::string_view keyName, DelphesParameters &parametersBlock) const
{
  std::ostringstream globalKeyName;
  globalKeyName << "::" << keyName;
  for(const std::string &childName : ChildrenList(globalKeyName.str()))
  {
    Tcl_Obj *variableName = Tcl_NewStringObj(const_cast<char *>(childName.data()), -1); // build the key name object
    Tcl_Obj *object = Tcl_ObjGetVar2(fTclInterp.get(), variableName, 0, TCL_GLOBAL_ONLY);
    if(const std::vector<std::string> subCollections = ChildrenList(childName); !subCollections.empty())
    {
      DelphesParameters subParameters;
      for(const std::string &subKeyName : subCollections)
        ParseParameters(subKeyName, subParameters);
      parametersBlock.Set(TrimmedName(childName), subParameters);
    }
    else
      ParseValue(object, childName, TrimmedName(childName), parametersBlock);
  }
}

void DelphesTCLConfReader::SetModuleType(std::string_view moduleName, std::string_view moduleType)
{
  if(const std::string moduleNameStr{moduleName}; fModuleTypes.count(moduleNameStr) > 0)
  {
    std::ostringstream message;
    message << "Duplicated module definition for module name '" << moduleName << "' with type '" << moduleType << "'.";
    throw std::runtime_error(message.str());
  }
  else
    fModuleTypes[moduleNameStr] = moduleType;
}

//------------------------------------------------------------------------------

int ModuleObjCmdProc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[])
{
  if(objc < 3)
  {
    Tcl_WrongNumArgs(interp, 1, objv, const_cast<char *>("className moduleName ?arg...?"));
    return TCL_ERROR;
  }

  // add module to a list of modules to be created
  DelphesTCLConfReader *reader = static_cast<DelphesTCLConfReader *>(clientData);
  reader->SetModuleType(Tcl_GetStringFromObj(objv[2], 0), Tcl_GetStringFromObj(objv[1], 0));

  if(objc > 3)
  {
    Tcl_Obj *object = Tcl_NewListObj(0, 0);
    Tcl_ListObjAppendElement(interp, object, Tcl_NewStringObj(const_cast<char *>("namespace"), -1));
    Tcl_ListObjAppendElement(interp, object, Tcl_NewStringObj(const_cast<char *>("eval"), -1));
    Tcl_ListObjAppendList(interp, object, Tcl_NewListObj(objc - 2, objv + 2));
    return Tcl_GlobalEvalObj(interp, object);
  }
  return TCL_OK;
}

//------------------------------------------------------------------------------

int SourceObjCmdProc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[])
{
  DelphesTCLConfReader *reader = static_cast<DelphesTCLConfReader *>(clientData);
  stringstream fileName;

  if(objc != 2)
  {
    Tcl_WrongNumArgs(interp, 1, objv, const_cast<char *>("fileName"));
    return TCL_ERROR;
  }

  fileName << Tcl_GetStringFromObj(objv[1], 0); //TODO: check if the base path is correct
  reader->ReadFile(fileName.str().data());

  return TCL_OK;
}

//------------------------------------------------------------------------------
