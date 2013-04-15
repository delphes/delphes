#ifndef ExRootTask_h
#define ExRootTask_h

/** \class ExRootTask
 *
 *  Class handling output ROOT tree
 *
 *  $Date: 2008-06-04 13:57:26 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "TTask.h"

#include "ExRootAnalysis/ExRootConfReader.h"

class TClass;
class TFolder;

class ExRootTask : public TTask
{
public:

  ExRootTask();
  virtual ~ExRootTask();

  virtual void Init();
  virtual void Process();
  virtual void Finish();

  virtual void InitTask();
  virtual void ProcessTask();
  virtual void FinishTask();

  virtual void InitSubTasks();
  virtual void ProcessSubTasks();
  virtual void FinishSubTasks();

  void Add(TTask *task);

  ExRootTask *NewTask(TClass *cl, const char *name);
  ExRootTask *NewTask(const char *className, const char *taskName);

  void Exec(Option_t* option);

  int GetInt(const char *name, int defaultValue, int index = -1);
  long GetLong(const char *name, long defaultValue, int index = -1);
  double GetDouble(const char *name, double defaultValue, int index = -1);
  bool GetBool(const char *name, bool defaultValue, int index = -1);
  const char *GetString(const char *name, const char *defaultValue, int index = -1);
  ExRootConfParam GetParam(const char *name);
  const ExRootConfReader::ExRootTaskMap *GetModules();

  void SetFolder(TFolder *folder) { fFolder = folder; }
  void SetConfReader(ExRootConfReader *conf) { fConfReader = conf; }

protected:

  TFolder *GetFolder() const { return fFolder; }
  ExRootConfReader *GetConfReader() const { return fConfReader; }

  TFolder *NewFolder(const char *name);
  TObject *GetObject(const char *name, TClass *cl);

private:

  TFolder *fFolder; //!
  ExRootConfReader *fConfReader; //!

  ClassDef(ExRootTask, 1)
};

#endif /* ExRootTask */

