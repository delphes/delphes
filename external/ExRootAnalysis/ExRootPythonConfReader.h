#ifndef ExRootPythonConfReader_h
#define ExRootPythonConfReader_h

/** \class ExRootPythonConfReader
 *
 *  Class handling input Python configuration steering
 *
 *  \author L. Forthomme - AGH, Kraków
 *
 */

#include "ExRootConfReader.h"

struct PyConfig;
struct _object;

class ExRootPythonConfParam : public ExRootConfParam
{
public:
  explicit ExRootPythonConfParam(const char *name = nullptr, _object *object = nullptr);

  int GetInt(int defaultValue) override;
  long GetLong(long defaultValue) override;
  double GetDouble(double defaultValue) override;
  bool GetBool(bool defaultValue) override;
  const char *GetString(const char *defaultValue) override;

  int GetSize() override;
  std::unique_ptr<ExRootConfParam> operator[](int index) override;

private:
  const char *fName{"\0"}; //!
  _object *fObject{nullptr};
};

//------------------------------------------------------------------------------

class ExRootPythonConfReader : public ExRootConfReader
{
public:
  ExRootPythonConfReader();
  ~ExRootPythonConfReader();

  void ReadFile(const char *fileName, bool isTop = true) override;
  std::unique_ptr<ExRootConfParam> GetParam(const char *name) override;

private:
  PyConfig *fConfig{nullptr};
  _object *fModule{nullptr};
  ClassDef(ExRootPythonConfReader, 1)
};

#endif
