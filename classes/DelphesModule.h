#ifndef DelphesModule_h
#define DelphesModule_h

/** \class DelphesModule
 *
 *  Base class for all Delphes modules
 *
 *  $Date: 2008-06-04 13:57:25 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootTask.h"

class TClass;
class TObject;
class TFolder;
class TClonesArray;

class ExRootResult;
class ExRootTreeBranch;
class ExRootTreeWriter;

class DelphesFactory;

class DelphesModule: public ExRootTask 
{
public:

  DelphesModule();
  ~DelphesModule();

  virtual void Init();
  virtual void Process();
  virtual void Finish();

  TObjArray *ImportArray(const char *name);
  TObjArray *ExportArray(const char *name);

  ExRootTreeBranch *NewBranch(const char *name, TClass *cl);

  ExRootResult *GetPlots();
  DelphesFactory *GetFactory();

protected:

  ExRootTreeWriter *fTreeWriter;
  DelphesFactory *fFactory;

private:

  ExRootResult *fPlots;

  TFolder *fPlotFolder, *fExportFolder;

  ClassDef(DelphesModule, 1)
};

#endif /* DelphesModule_h */

