#ifndef ExRootConfReader_h
#define ExRootConfReader_h

/** \class ExRootConfReader
 *
 *  Class handling configuration data
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "TNamed.h"

struct Jim_Obj;
struct Jim_Interp;

class ExRootConfParam: public TNamed
{
public:
  ExRootConfParam(const char *name = 0, Jim_Obj *object = 0, Jim_Interp *interp = 0);

  int GetInt(int defaultValue = 0);
  long GetLong(long defaultValue = 0);
  double GetDouble(double defaultValue = 0.0);
  bool GetBool(bool defaultValue = false);
  const char *GetString(const char *defaultValue = "");

  int GetSize();
  ExRootConfParam operator[](int index);

private:
  Jim_Obj *fObject; //!
  Jim_Interp *fTclInterp; //!

  ClassDef(ExRootConfParam, 1)
};

//------------------------------------------------------------------------------

class ExRootConfReader: public TNamed
{
public:
  ExRootConfReader();
  ~ExRootConfReader();

  void ReadData(const char *dirName, char *data, int length);
  void ReadFile(const char *fileName, bool isTop = true);

  int GetInt(const char *name, int defaultValue, int index = -1);
  long GetLong(const char *name, long defaultValue, int index = -1);
  double GetDouble(const char *name, double defaultValue, int index = -1);
  bool GetBool(const char *name, bool defaultValue, int index = -1);
  const char *GetString(const char *name, const char *defaultValue, int index = -1);
  ExRootConfParam GetParam(const char *name);

  const char *GetTopDir() const { return fTopDir; }

private:
  const char *fTopDir; //!

  Jim_Interp *fTclInterp; //!

  ClassDef(ExRootConfReader, 1)
};

#endif
