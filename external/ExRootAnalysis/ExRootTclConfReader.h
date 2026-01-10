#ifndef ExRootTclConfReader_h
#define ExRootTclConfReader_h

/** \class ExRootTclConfReader
 *
 *  Class handling input Tcl configuration steering
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootConfReader.h"

struct Tcl_Obj;
struct Tcl_Interp;

class ExRootTclConfParam : public ExRootConfParam
{
public:
  explicit ExRootTclConfParam(const char *name = 0, Tcl_Obj *object = 0, Tcl_Interp *interp = 0);

  int GetInt(int defaultValue) override;
  long GetLong(long defaultValue) override;
  double GetDouble(double defaultValue) override;
  bool GetBool(bool defaultValue) override;
  const char *GetString(const char *defaultValue) override;

  int GetSize() override;
  std::unique_ptr<ExRootConfParam> operator[](int index) override;

private:
  const char *fName; //!
  Tcl_Obj *fObject; //!
  Tcl_Interp *fTclInterp; //!
};

//------------------------------------------------------------------------------

class ExRootTclConfReader : public ExRootConfReader
{
public:
  ExRootTclConfReader();
  ~ExRootTclConfReader();

  void ReadFile(const char *fileName, bool isTop = true) override;
  std::unique_ptr<ExRootConfParam> GetParam(const char *name) override;
  std::unique_ptr<ExRootConfParam> GetGlobalParam(const char *name) override;

private:
  Tcl_Interp *fTclInterp; //!

  ClassDef(ExRootTclConfReader, 1)
};

#endif
