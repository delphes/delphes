#ifndef ExRootConfReader_h
#define ExRootConfReader_h

/** \class ExRootConfReader
 *
 *  Class handling output ROOT tree
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "TNamed.h"

#include <map>
#include <utility>

struct Tcl_Obj;
struct Tcl_Interp;

class ExRootConfParam
{
public:

  ExRootConfParam(const char *name = 0, Tcl_Obj *object = 0, Tcl_Interp *interp = 0);

  int GetInt(int defaultValue = 0);
  long GetLong(long defaultValue = 0);
  double GetDouble(double defaultValue = 0.0);
  bool GetBool(bool defaultValue = false);
  const char *GetString(const char *defaultValue = "");

  int GetSize();
  ExRootConfParam operator[](int index);

private:

  const char *fName; //!
  Tcl_Obj *fObject; //!
  Tcl_Interp *fTclInterp; //!
};

//------------------------------------------------------------------------------

class ExRootConfReader : public TNamed
{
public:
  typedef std::map<TString, TString> ExRootTaskMap;

  ExRootConfReader();
  ~ExRootConfReader();

  void ReadFile(const char *fileName);

  int GetInt(const char *name, int defaultValue, int index = -1);
  long GetLong(const char *name, long defaultValue, int index = -1);
  double GetDouble(const char *name, double defaultValue, int index = -1);
  bool GetBool(const char *name, bool defaultValue, int index = -1);
  const char *GetString(const char *name, const char *defaultValue, int index = -1);
  ExRootConfParam GetParam(const char *name);

  const ExRootTaskMap *GetModules() const { return &fModules; }

  void AddModule(const char *className, const char *moduleName);

private:

  Tcl_Interp *fTclInterp; //!

  ExRootTaskMap fModules; //!

  ClassDef(ExRootConfReader, 1)
};

#endif

