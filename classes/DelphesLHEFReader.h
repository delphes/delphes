#ifndef DelphesLHEFReader_h
#define DelphesLHEFReader_h

/** \class DelphesLHEFReader
 *
 *  Reads LHEF file
 *
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include <stdio.h>

#include <vector>

class TObjArray;
class TStopwatch;
class TDatabasePDG;
class ExRootTreeBranch;
class DelphesFactory;

class DelphesLHEFReader
{
public:

  DelphesLHEFReader();
  ~DelphesLHEFReader();

  void SetInputFile(FILE *inputFile);

  void Clear();
  bool EventReady();

  bool ReadBlock(DelphesFactory *factory,
    TObjArray *allParticleOutputArray,
    TObjArray *stableParticleOutputArray,
    TObjArray *partonOutputArray);

  void AnalyzeEvent(ExRootTreeBranch *branch, long long eventNumber,
    TStopwatch *readStopWatch, TStopwatch *procStopWatch);

  void AnalyzeRwgt(ExRootTreeBranch *branch);

private:

  void AnalyzeParticle(DelphesFactory *factory,
    TObjArray *allParticleOutputArray,
    TObjArray *stableParticleOutputArray,
    TObjArray *partonOutputArray);

  FILE *fInputFile;

  char *fBuffer;

  TDatabasePDG *fPDG;

  bool fEventReady;

  int fEventCounter;

  int fParticleCounter, fProcessID;
  double fWeight, fScalePDF, fAlphaQCD, fAlphaQED;

  int fPID, fStatus, fM1, fM2, fC1, fC2;
  double fPx, fPy, fPz, fE, fMass;
  
  std::vector<double> fRwgtList;
};

#endif // DelphesLHEFReader_h


