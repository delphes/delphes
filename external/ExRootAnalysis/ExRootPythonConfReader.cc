
/** \class ExRootConfReader
 *
 *  Class handling output ROOT tree
 *
 *  \author L. Forthomme - AGH, Kraków
 *
 */

#include "ExRootAnalysis/ExRootPythonConfReader.h"

#include <Python.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TSystem.h>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

//------------------------------------------------------------------------------

ExRootPythonConfReader::ExRootPythonConfReader() : fConfig(new PyConfig)
{
}

//------------------------------------------------------------------------------

ExRootPythonConfReader::~ExRootPythonConfReader()
{
  if(fModule) Py_DECREF(fModule);
  if(fConfig) delete fConfig;
}

//------------------------------------------------------------------------------

void ExRootPythonConfReader::ReadFile(const char *fileName, bool isTop)
{
  if(isTop) fTopDir = gSystem->DirName(fileName);

#ifdef _WIN32
  static constexpr char PATH_DELIM = ';';
#else
  static constexpr char PATH_DELIM = ':';
#endif

  std::ostringstream os;
  os << fTopDir;
  if(const auto *delphes_path = getenv("Delphes_PATH"); delphes_path)
    os << PATH_DELIM << delphes_path;
  os << PATH_DELIM << "../python";
  setenv("PYTHONPATH", os.str().data(), 1);
  std::cout << getenv("PYTHONPATH") << std::endl;
  const auto file_path = std::string{gSystem->BaseName(fileName)};
#if PY_VERSION_HEX >= 0x03080000
  PyConfig_InitPythonConfig(fConfig);
  fConfig->parser_debug = 0;
  fConfig->verbose = 0;
  Py_InitializeFromConfig(fConfig);
#else
  Py_DebugFlag = 0;
  Py_VerboseFlag = 0;
  Py_InitializeEx(1);
#endif
  const auto rawname = file_path.substr(0, file_path.find_last_of("."));
  if(fModule = PyImport_ImportModule(rawname.c_str()); !fModule)
  {
    std::ostringstream message;
    message << "Failed to read Python configuration file at '" << fileName << "'" << std::endl;
    PyErr_Print();
    throw std::runtime_error(message.str());
  }
  const auto execution_path = GetGlobalParam("ExecutionPath");
  if(!PyObject_HasAttrString(fModule, "__dir__"))
  {
    std::ostringstream message;
    message << "Failed to retrieve the Python attributes search method for the configuration file at '" << fileName << "'" << std::endl;
    throw std::runtime_error(message.str());
  }
  for(long i = 0; i < execution_path->GetSize(); ++i)
  {
    const auto module_name = (*execution_path)[i]->GetString(),
               module_type = GetGlobalParam(module_name)->GetParam("type")->GetString();
    AddModule(module_type, module_name); //TODO: retrieve it from parameters name
  }
}

//------------------------------------------------------------------------------

std::unique_ptr<ExRootConfParam> ExRootPythonConfReader::GetParam(const char *name)
{
  const auto path = TString(name).Tokenize("::");
  if(path->GetEntriesFast() <= 1)
  {
    if(PyObject_HasAttrString(fModule, name) == 1)
      return std::make_unique<ExRootPythonConfParam>(name, PyObject_GetAttrString(fModule, name), true);
    std::cout << "WARNING: Failed to retrieve a parameter with name '" << name << "'." << std::endl;
    return std::make_unique<ExRootPythonConfParam>(name, nullptr, true);
  }
  else
  {
    std::unique_ptr<ExRootConfParam> out;
    for(int j = 0; j < path->GetEntriesFast(); j++)
    {
      const auto dir = reinterpret_cast<TObjString *>(path->At(j))->GetString();
      if(j == 0)
        out = GetParam(dir);
      else
        out = out->GetParam(dir);
    }
    return out;
  }
  std::cout << "WARNING: Failed to retrieve a parameter with name '" << name << "'." << std::endl;
  return std::make_unique<ExRootPythonConfParam>(name, nullptr, true);
}

//------------------------------------------------------------------------------

ExRootPythonConfParam::ExRootPythonConfParam(const char *name, PyObject *object, bool wrapOnly) :
  ExRootConfParam(name), fObject(object), fWrapOnly(wrapOnly)
{
}

//------------------------------------------------------------------------------

ExRootPythonConfParam::~ExRootPythonConfParam()
{
  if(fObject && !fWrapOnly) Py_DECREF(fObject);
}

//------------------------------------------------------------------------------

int ExRootPythonConfParam::GetInt(int defaultValue)
{
  if(!fObject)
    return defaultValue;
  if(fObject && !PyLong_Check(fObject))
  {
    std::ostringstream message;
    message << "Parameter '" << fName << "' is not a Python integer." << std::endl;
    throw std::runtime_error(message.str());
  }
  return PyLong_AsLongLong(fObject);
}

//------------------------------------------------------------------------------

long ExRootPythonConfParam::GetLong(long defaultValue)
{
  if(!fObject)
    return defaultValue;
  if(fObject && !PyLong_Check(fObject))
  {
    std::ostringstream message;
    message << "Parameter '" << fName << "' is not a Python long integer." << std::endl;
    throw std::runtime_error(message.str());
  }
  return PyLong_AsLongLong(fObject);
}

//------------------------------------------------------------------------------

double ExRootPythonConfParam::GetDouble(double defaultValue)
{
  if(!fObject)
    return defaultValue;
  if(!PyFloat_Check(fObject))
  {
    std::ostringstream message;
    message << "parameter '" << fName << "' is not a Python float." << std::endl;
    throw std::runtime_error(message.str());
  }
  return PyFloat_AsDouble(fObject);
}

//------------------------------------------------------------------------------

bool ExRootPythonConfParam::GetBool(bool defaultValue)
{
  if(!fObject)
    return defaultValue;
  if(!PyBool_Check(fObject))
  {
    std::ostringstream message;
    message << "parameter '" << fName << "' is not a Python boolean." << std::endl;
    throw std::runtime_error(message.str());
  }
  return PyObject_IsTrue(fObject);
}

//------------------------------------------------------------------------------

const char *ExRootPythonConfParam::GetString(const char *defaultValue)
{
  if(!fObject)
    return defaultValue;
  if(PyUnicode_Check(fObject))
    return PyUnicode_AsUTF8(fObject);
  if(PyBytes_Check(fObject))
    return ::strdup(PyBytes_AS_STRING(fObject));
  std::ostringstream message;
  message << "parameter '" << fName << "' is not a Python string." << std::endl;
  throw std::runtime_error(message.str());
}

//------------------------------------------------------------------------------

std::unique_ptr<ExRootConfParam> ExRootPythonConfParam::GetParam(const char *paramName) const
{
  if(!fObject)
  {
    std::ostringstream message;
    message << "Parameter '" << fName << "' is not a valid object. "
            << "Cannot retrieve a sub-parameters '" << paramName << "' from it." << std::endl;
    throw std::runtime_error(message.str());
  }
  if(PyDict_Check(fObject) == 1) // we are in a dictionary, we can try to retrieve an element with key
    return std::make_unique<ExRootPythonConfParam>(paramName,
      PyDict_GetItemString(fObject, paramName),
      true);
  if(PyObject_HasAttrString(fObject, paramName) == 1) // retrieve an element from global scope
    return std::make_unique<ExRootPythonConfParam>(paramName,
      PyObject_GetAttrString(fObject, paramName),
      true);
  std::cout << "WARNING: Failed to retrieve a parameter with name '" << paramName << "' from object "
            << "'" << fName << "'." << std::endl;
  return std::make_unique<ExRootPythonConfParam>(paramName, nullptr, true);
}

//------------------------------------------------------------------------------

int ExRootPythonConfParam::GetSize()
{
  if(!fObject)
  {
    std::cout << "WARNING: Parameter '" << fName << "' is not a valid object. "
              << "Cannot retrieve its size." << std::endl;
    return 0;
  }
  if(PySequence_Check(fObject))
    return PySequence_Size(fObject);
  std::ostringstream message;
  message << "Cannot retrieve object size. Parameter '" << fName << "' is not a sequence." << std::endl;
  throw std::runtime_error(message.str());
}

//------------------------------------------------------------------------------

std::unique_ptr<ExRootConfParam> ExRootPythonConfParam::operator[](int index)
{
  if(!fObject)
  {
    std::ostringstream message;
    message << "parameter '" << fName << "' is not a valid object." << std::endl;
    throw std::runtime_error(message.str());
  }
  if(PySequence_Check(fObject))
    return std::make_unique<ExRootPythonConfParam>(fName, PySequence_GetItem(fObject, index), true);
  std::ostringstream message;
  message << "parameter '" << fName << "' is not a sequence." << std::endl;
  throw std::runtime_error(message.str());
}
