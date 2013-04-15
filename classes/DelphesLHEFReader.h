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

private:

  void AnalyzeParticle(DelphesFactory *factory,
    TObjArray *allParticleOutputArray,
    TObjArray *stableParticleOutputArray,
    TObjArray *partonOutputArray);

  FILE *fInputFile;

  char *fBuffer;

  TDatabasePDG *fPDG;

  int fEventCounter;

  int fParticleCounter, fProcessID;
  double fWeight, fScalePDF, fAlphaQCD, fAlphaQED;

  int fPID, fStatus, fM1, fM2, fC1, fC2;
  double fPx, fPy, fPz, fE, fMass;
};

#endif // DelphesLHEFReader_h


