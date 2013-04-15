#ifndef DelphesHepMCReader_h
#define DelphesHepMCReader_h

/** \class DelphesHepMCReader
 *
 *  Reads HepMC file
 *
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include <map>

#include <stdio.h>

class TObjArray;
class TStopwatch;
class TDatabasePDG;
class ExRootTreeBranch;
class DelphesFactory;

class DelphesHepMCReader
{
public:

  DelphesHepMCReader();
  ~DelphesHepMCReader();

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

  void FinalizeParticles(TObjArray *allParticleOutputArray);

  FILE *fInputFile;

  char *fBuffer;

  TDatabasePDG *fPDG;

  int fEventNumber, fMPI, fProcessID, fSignalCode, fVertexCounter;
  double fScale, fAlphaQCD, fAlphaQED;

  int fID1, fID2;
  double fX1, fX2, fScalePDF, fPDF1, fPDF2;

  int fOutVertexCode, fVertexID, fInCounter, fOutCounter;
  double fX, fY, fZ, fT;

  int fParticleCode, fPID, fStatus, fInVertexCode;
  double fPx, fPy, fPz, fE, fMass, fTheta, fPhi;

  int fParticleCounter;

  std::map< int, std::pair < int, int > > fMotherMap;
  std::map< int, std::pair < int, int > > fDaughterMap;
};

#endif // DelphesHepMCReader_h


