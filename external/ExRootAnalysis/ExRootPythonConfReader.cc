
/** \class ExRootConfReader
 *
 *  Class handling output ROOT tree
 *
 *  \author L. Forthomme - AGH, Kraków
 *
 */

#include "ExRootAnalysis/ExRootPythonConfReader.h"

#include <Python.h>
#include <TSystem.h>

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

//------------------------------------------------------------------------------

ExRootPythonConfReader::ExRootPythonConfReader() : fConfig(new PyConfig)
{
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

  if(fModule = PyImport_AddModule(fileName); !fModule)
  {
    std::ostringstream message;
    message << "can't read configuration file " << fileName << std::endl;
    throw std::runtime_error(message.str());
  }
}

//------------------------------------------------------------------------------

std::unique_ptr<ExRootConfParam> ExRootPythonConfReader::GetParam(const char *name)
{
  if(PyObject_HasAttrString(fModule, name) == 1)
    return std::make_unique<ExRootPythonConfParam>(name, PyObject_GetAttrString(fModule, name));
  return std::make_unique<ExRootPythonConfParam>(name, nullptr);
}

//------------------------------------------------------------------------------

ExRootPythonConfParam::ExRootPythonConfParam(const char *name, PyObject *object) :
  ExRootConfParam(name), fObject(object)
{
}

//------------------------------------------------------------------------------

int ExRootPythonConfParam::GetInt(int defaultValue)
{
  if(!fObject)
    return defaultValue;
  if(fObject && !PyLong_Check(fObject))
  {
    std::ostringstream message;
    message << "parameter '" << fName << "' is not an integer." << std::endl;
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
    message << "parameter '" << fName << "' is not an long integer." << std::endl;
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
    message << "parameter '" << fName << "' is not a number." << std::endl;
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
    message << "parameter '" << fName << "' is not a boolean." << std::endl;
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
    if(auto *str_buf = ::strdup(PyBytes_AS_STRING(fObject)); str_buf)
    {
      auto out = std::string(str_buf);
      free(str_buf);
      return out.data();
    }
  std::ostringstream message;
  message << "parameter '" << fName << "' is not a string." << std::endl;
  throw std::runtime_error(message.str());
}

//------------------------------------------------------------------------------

int ExRootPythonConfParam::GetSize()
{
  if(!fObject)
  {
    std::ostringstream message;
    message << "parameter '" << fName << "' is not a valid object." << std::endl;
    throw std::runtime_error(message.str());
  }
  if(PyTuple_Check(fObject))
    return PyTuple_Size(fObject);
  if(PyList_Check(fObject))
    return PyList_Size(fObject);
  std::ostringstream message;
  message << "parameter '" << fName << "' is not a tuple nor a list." << std::endl;
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
  if(PyTuple_Check(fObject))
    return std::make_unique<ExRootPythonConfParam>(fName, PyTuple_GetItem(fObject, index));
  if(PyList_Check(fObject))
    return std::make_unique<ExRootPythonConfParam>(fName, PyList_GetItem(fObject, index));
  std::ostringstream message;
  message << "parameter '" << fName << "' is not a tuple nor a list." << std::endl;
  throw std::runtime_error(message.str());
}
