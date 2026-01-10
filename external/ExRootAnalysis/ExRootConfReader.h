#ifndef ExRootConfReader_h
#define ExRootConfReader_h

/** \class ExRootConfReader
 *
 *  Class handling input steering card
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "TNamed.h"

#include <map>
#include <memory>
#include <utility>

class ExRootConfParam
{
public:
  explicit ExRootConfParam(const char *name = 0);
  virtual ~ExRootConfParam() = default;

  virtual int GetInt(int defaultValue = 0) = 0;
  virtual long GetLong(long defaultValue = 0) = 0;
  virtual double GetDouble(double defaultValue = 0.0) = 0;
  virtual bool GetBool(bool defaultValue = false) = 0;
  virtual const char *GetString(const char *defaultValue = "") = 0;

  virtual int GetSize() = 0;
  virtual std::unique_ptr<ExRootConfParam> operator[](int index) = 0;

protected:
  const char *fName; //!
};

//------------------------------------------------------------------------------

class ExRootConfReader : public TNamed
{
public:
  typedef std::map<TString, TString> ExRootTaskMap;

  static std::unique_ptr<ExRootConfReader> ReadConf(const char *filename);

  virtual ~ExRootConfReader() = default;

  virtual void ReadFile(const char *fileName, bool isTop = true) = 0;

  int GetInt(const char *name, int defaultValue, int index = -1);
  long GetLong(const char *name, long defaultValue, int index = -1);
  double GetDouble(const char *name, double defaultValue, int index = -1);
  bool GetBool(const char *name, bool defaultValue, int index = -1);
  const char *GetString(const char *name, const char *defaultValue, int index = -1);

  virtual std::unique_ptr<ExRootConfParam> GetParam(const char *name) = 0;
  virtual std::unique_ptr<ExRootConfParam> GetGlobalParam(const char *name);

  const ExRootTaskMap *GetModules() const { return &fModules; }

  void AddModule(const char *className, const char *moduleName);

  const char *GetTopDir() const { return fTopDir; }

protected:
  ExRootConfReader();

  const char *fTopDir; //!
  ExRootTaskMap fModules; //!
};

#endif
