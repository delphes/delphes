#ifndef DelphesSTDHEPReader_h
#define DelphesSTDHEPReader_h

/** \class DelphesSTDHEPReader
 *
 *  Reads STDHEP file
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
#include <rpc/types.h>
#include <rpc/xdr.h>

class TObjArray;
class TStopwatch;
class TDatabasePDG;
class ExRootTreeBranch;
class DelphesFactory;

class DelphesSTDHEPReader
{
public:
  enum STDHEPBlock
  {
    GENERIC = 0,
    FILEHEADER = 1,
    EVENTTABLE = 2,
    EVENTHEADER = 4,
    MCFIO_STDHEP = 101,
    MCFIO_STDHEPBEG = 106,
    MCFIO_STDHEPEND = 107,
    MCFIO_STDHEP4 = 201
  };

  DelphesSTDHEPReader();
  ~DelphesSTDHEPReader();

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

  void AnalyzeParticles(DelphesFactory *factory,
    TObjArray *allParticleOutputArray,
    TObjArray *stableParticleOutputArray,
    TObjArray *partonOutputArray);

  void SkipBytes(u_int size);
  void SkipArray(u_int elsize);

  void ReadFileHeader();
  void ReadEventTable();
  void ReadEventHeader();
  void ReadSTDCM1();
  void ReadSTDHEP();
  void ReadSTDHEP4();

  FILE *fInputFile;

  XDR *fInputXDR;

  char *fBuffer;

  TDatabasePDG *fPDG;

  u_int fEntries;
  int fBlockType, fEventNumber, fEventSize;
  double fWeight, fAlphaQCD, fAlphaQED;

  u_int fScaleSize;
  double fScale[10];
};

#endif // DelphesSTDHEPReader_h


