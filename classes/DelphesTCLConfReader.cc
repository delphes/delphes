/** \class DelphesTCLConfReader
 *
 *  Class handling TCL input card parsing
 *
 *  \author L. Forthomme - AGH, Krakow
 *  \note Based on ExRootConfReader by P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesTCLConfReader.h"

#include <tcl/tcl.h>

#include <fstream>

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

void DelphesTCLConfReader::ReadFile(std::string_view fileName)
{
  std::ifstream inputFileStream(std::string{fileName}, std::ios::in);
  if(!inputFileStream.is_open())
  {
    std::ostringstream message;
    message << "can't open configuration file " << fileName;
    throw std::runtime_error(message.str());
  }
  std::string buffer((std::istreambuf_iterator<char>(inputFileStream)), std::istreambuf_iterator<char>());

  Tcl_Obj *cmdObjPtr = Tcl_NewObj();
  cmdObjPtr->bytes = buffer.data();
  cmdObjPtr->length = buffer.size();

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
  ParseValue(nullptr, "::MaxEvents", "MaxEvents", fParams);
  ParseValue(nullptr, "::RandomSeed", "RandomSeed", fParams);
  ParseValue(nullptr, "::SkipEvents", "SkipEvents", fParams);

  for(const std::string &moduleName : fParams.Get<std::vector<std::string> >("ExecutionPath"))
  {
    DelphesParameters moduleParams;
    ParseParameters(moduleName, moduleParams);
    moduleParams.Set("ModuleType", fModuleTypes.at(moduleName));
    fParams.Set(moduleName, moduleParams);
  }

  std::ofstream outputCard("lastrun.yaml");
  outputCard << fParams;
}

//------------------------------------------------------------------------------

template <>
std::string DelphesTCLConfReader::Get<std::string>(Tcl_Obj *objPtr) const
{
  return Tcl_GetStringFromObj(objPtr, 0);
}

//------------------------------------------------------------------------------

template <>
int DelphesTCLConfReader::Get<int>(Tcl_Obj *objPtr) const
{
  if(int result; Tcl_GetIntFromObj(fTclInterp.get(), objPtr, &result) == TCL_OK)
    return result;
  std::ostringstream message;
  message << "Failed to cast object at " << objPtr << " to integer.";
  if(Tcl_Obj *error = Tcl_GetObjResult(fTclInterp.get()); error)
    message << " Tcl error: " << Get<std::string>(error);
  throw std::runtime_error(message.str());
}

//------------------------------------------------------------------------------

template <>
double DelphesTCLConfReader::Get<double>(Tcl_Obj *objPtr) const
{
  if(double result; Tcl_GetDoubleFromObj(fTclInterp.get(), objPtr, &result) == TCL_OK)
    return result;
  std::ostringstream message;
  message << "Failed to cast object at " << objPtr << " to double floating point number.";
  if(Tcl_Obj *error = Tcl_GetObjResult(fTclInterp.get()); error)
    message << " Tcl error: " << Get<std::string>(error);
  throw std::runtime_error(message.str());
}

//------------------------------------------------------------------------------

std::vector<Tcl_Obj *> DelphesTCLConfReader::GetObjVector(Tcl_Obj *objPtr) const
{
  Tcl_Obj **objVectorPtr = nullptr;
  if(int length; Tcl_ListObjGetElements(fTclInterp.get(), objPtr, &length, &objVectorPtr) == TCL_OK)
    return std::vector<Tcl_Obj *>(objVectorPtr, objVectorPtr + length);
  std::ostringstream message;
  message << "Failed to retrieve length of object at " << objPtr << ". Is it really a vector?";
  if(Tcl_Obj *error = Tcl_GetObjResult(fTclInterp.get()); error)
    message << " Error: " << Get<std::string>(error);
  throw std::runtime_error(message.str());
}

//------------------------------------------------------------------------------

std::string DelphesTCLConfReader::Run(std::string_view command) const
{
  if(Tcl_Eval(fTclInterp.get(), const_cast<char *>(command.data())) != TCL_OK) return "nil";
  return Get<std::string>(Tcl_GetObjResult(fTclInterp.get()));
}

//------------------------------------------------------------------------------

void DelphesTCLConfReader::ParseValue(const Tcl_Obj *tclObject,
  std::string_view tclName, std::string_view keyName,
  DelphesParameters &delphesParams) const
{
  Tcl_Obj *object = nullptr;
  if(tclObject)
    object = const_cast<Tcl_Obj *>(tclObject);
  else
  { // retrieve from the global environment
    if(Tcl_Obj *varNameObj = Tcl_NewStringObj(const_cast<char *>(tclName.data()), -1); varNameObj) // build the key name object
      object = Tcl_ObjGetVar2(fTclInterp.get(), varNameObj, 0, TCL_GLOBAL_ONLY);
    else
      throw std::runtime_error("Failed to build a Tcl string for the variable name '" + std::string{tclName} + "'. "
        + "Please check the status of your interpreter.");
  }
  if(!object)
  {
    //std::cout << "Failed to retrieve object with key '" << tclName << "' from object.";
    return;
  }

  const std::string delphesKeyName = keyName.empty() ? std::string{tclName} : std::string{keyName};

  if(const std::vector<Tcl_Obj *> collObj = GetObjVector(object); collObj.size() > 1)
  { // we have a list/collection
    // with this version of TCL, only flat lists are supported, no need to parse dictionaries
    std::vector<int> intCollection;
    std::vector<double> doubleCollection;
    std::vector<bool> boolCollection;
    std::vector<std::string> stringCollection;
    for(Tcl_Obj *const &subObject : collObj)
    {
      if(int result; Tcl_GetIntFromObj(fTclInterp.get(), subObject, &result) == TCL_OK)
        intCollection.emplace_back(result);
      else if(double result; Tcl_GetDoubleFromObj(fTclInterp.get(), subObject, &result) == TCL_OK)
        doubleCollection.emplace_back(result);
      else if(int result; Tcl_GetBooleanFromObj(fTclInterp.get(), object, &result) == TCL_OK)
        boolCollection.emplace_back(result);
      else if(const std::string result = Get<std::string>(subObject); !result.empty())
        stringCollection.emplace_back(result);
    }

    //------------------------------------------------------------------------------
    // here starts the hacky part to match TCL cards to structured configurations
    //------------------------------------------------------------------------------

    if(delphesKeyName == "Branch"
      && stringCollection.size() > 1 && stringCollection.size() % 3 == 0)
    {
      std::vector<std::array<std::string, 3> > branchesInfo;
      for(size_t i = 0; i < stringCollection.size() / 3; ++i)
        branchesInfo.emplace_back(std::array{
          stringCollection.at(3 * i), stringCollection.at(3 * i + 1), stringCollection.at(3 * i + 2)});
      delphesParams.Set(delphesKeyName, branchesInfo);
      return;
    }
    if(delphesKeyName == "EnergyFraction")
    {
      std::unordered_map<unsigned long long, double> singleEnergyFraction;
      std::unordered_map<unsigned long long, std::pair<double, double> > doubleEnergyFraction;
      const std::vector<Tcl_Obj *> objColl = GetObjVector(object);
      // read energy fractions for different particles
      for(size_t i = 0; i < objColl.size() / 2; ++i)
      {
        const int pdgId = Get<int>(objColl.at(2 * i));
        if(const std::vector<double> energyFractions = GetVector<double>(objColl.at(2 * i + 1)); energyFractions.size() == 1)
          singleEnergyFraction[pdgId] = energyFractions.at(0);
        else if(energyFractions.size() == 2)
          doubleEnergyFraction[pdgId] = std::make_pair(energyFractions.at(0), energyFractions.at(1));
      }
      if(!singleEnergyFraction.empty())
        delphesParams.Set(delphesKeyName, singleEnergyFraction);
      else if(!doubleEnergyFraction.empty())
        delphesParams.Set(delphesKeyName, doubleEnergyFraction);
      return;
    }
    if(delphesKeyName.find("EtaPhiBins") != std::string::npos)
    {
      std::unordered_map<double, std::vector<double> > etaPhiBins;
      const std::vector<Tcl_Obj *> objColl = GetObjVector(object);
      for(size_t i = 0; i < objColl.size() / 2; ++i)
      {
        const std::vector<Tcl_Obj *> phiBinsObj = GetObjVector(objColl.at(2 * i + 1)); // phi binning
        for(Tcl_Obj *const &etaValueObj : GetObjVector(objColl.at(2 * i))) // eta binning
        {
          const double etaValue = Get<double>(etaValueObj);
          for(Tcl_Obj *const &phiValueObj : phiBinsObj)
            etaPhiBins[etaValue].emplace_back(Get<double>(phiValueObj));
        }
      }
      delphesParams.Set(delphesKeyName, etaPhiBins);
      return;
    }
    if(delphesKeyName.find("Formula") != std::string::npos)
    {
      if(!intCollection.empty() && doubleCollection.size() == intCollection.size())
      { //FIXME: hacky way to ensure {pdgId -> formula} are parsed correctly
        std::unordered_map<long, std::string> intBasedFormula;
        for(size_t i = 0; i < intCollection.size(); ++i)
          intBasedFormula[intCollection.at(i)] = std::to_string(doubleCollection.at(i));
        delphesParams.Set(delphesKeyName, intBasedFormula);
        return;
      }
      if(!intCollection.empty() && stringCollection.size() == intCollection.size())
      { //FIXME: hacky way to ensure {pdgId -> formula} are parsed correctly
        std::unordered_map<long, std::string> intBasedFormula;
        for(size_t i = 0; i < intCollection.size(); ++i)
          intBasedFormula[intCollection.at(i)] = stringCollection.at(i);
        delphesParams.Set(delphesKeyName, intBasedFormula);
        return;
      }
      { //FIXME: hacky way to ensure formulas are parsed as full strings
        if(const std::string result = Get<std::string>(object); !result.empty())
        { //FIXME check if reasonable
          delphesParams.Set(delphesKeyName, result);
          return;
        }
      }
    }
    if(delphesKeyName == "RhoEtaRange" // tweak for FastJetFinder module
      && !doubleCollection.empty() && doubleCollection.size() % 2 == 0)
    {
      std::vector<std::pair<double, double> > valuesRanges;
      for(size_t i = 0; i < doubleCollection.size() / 2; ++i)
        valuesRanges.emplace_back(std::make_pair(doubleCollection.at(2 * i), doubleCollection.at(2 * i + 1)));
      delphesParams.Set(delphesKeyName, valuesRanges);
      return;
    }
    if(delphesKeyName == "InputArray"
      && stringCollection.size() > 1 && stringCollection.size() % 2 == 0
      && stringCollection.at(1).find("/") == std::string::npos)
    {
      std::vector<std::pair<std::string, std::string> > inputOutputArrays;
      for(size_t i = 0; i < stringCollection.size() / 2; ++i)
        inputOutputArrays.emplace_back(std::make_pair(stringCollection.at(2 * i), stringCollection.at(2 * i + 1)));
      delphesParams.Set(delphesKeyName, inputOutputArrays);
      return;
    }
    if(delphesKeyName.find("Array") != std::string::npos && stringCollection.size() == 1)
    {
      delphesParams.Set(delphesKeyName, stringCollection.at(0));
      return;
    }

    //------------------------------------------------------------------------------
    // end of the hacky part, you may now breathe normally
    //------------------------------------------------------------------------------

    if(!intCollection.empty())
    {
      if(intCollection.size() == 1)
        delphesParams.Set(delphesKeyName, intCollection);
      else
        delphesParams.Set(delphesKeyName, intCollection);
    }
    if(!doubleCollection.empty())
    {
      if(doubleCollection.size() == 1)
        delphesParams.Set(delphesKeyName, doubleCollection.at(0));
      else
        delphesParams.Set(delphesKeyName, doubleCollection);
    }
    if(!boolCollection.empty())
    {
      if(boolCollection.size() == 1)
        delphesParams.Set(delphesKeyName, boolCollection);
      else
        delphesParams.Set(delphesKeyName, boolCollection);
    }
    if(!stringCollection.empty())
    {
      if(stringCollection.size() == 1)
        delphesParams.Set(delphesKeyName, stringCollection);
      else
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
  else if(const std::string result = Get<std::string>(object); !result.empty())
    delphesParams.Set(delphesKeyName, result); // fallback solution: parse everything as a string
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
      childrenList[i] = Get<std::string>(subObject);
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
      return Get<std::string>(trimmedSubObject);
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
  std::ostringstream fileName;

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
