#ifndef TreeWriter_h
#define TreeWriter_h

/** \class TreeWriter
 *
 *  Fills ROOT tree branches.
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <map>

class TClass;
class TObjArray;
class TRefArray;

class Candidate;
class ExRootTreeBranch;

class TreeWriter: public DelphesModule
{
public:

  TreeWriter();
  ~TreeWriter();

  void Init();
  void Process();
  void Finish();

private:

  void FillParticles(Candidate *candidate, TRefArray *array);

  void ProcessParticles(ExRootTreeBranch *branch, TObjArray *array);
  void ProcessVertices(ExRootTreeBranch *branch, TObjArray *array);
  void ProcessTracks(ExRootTreeBranch *branch, TObjArray *array);
  void ProcessTowers(ExRootTreeBranch *branch, TObjArray *array);
  void ProcessPhotons(ExRootTreeBranch *branch, TObjArray *array);
  void ProcessElectrons(ExRootTreeBranch *branch, TObjArray *array);
  void ProcessMuons(ExRootTreeBranch *branch, TObjArray *array);
  void ProcessTauJets(ExRootTreeBranch *branch, TObjArray *array);
  void ProcessJets(ExRootTreeBranch *branch, TObjArray *array);
  void ProcessMissingET(ExRootTreeBranch *branch, TObjArray *array);
  void ProcessScalarHT(ExRootTreeBranch *branch, TObjArray *array);
  void ProcessRho(ExRootTreeBranch *branch, TObjArray *array);
  void ProcessWeight(ExRootTreeBranch *branch, TObjArray *array);
  void ProcessHectorHit(ExRootTreeBranch *branch, TObjArray *array);

#if !defined(__CINT__) && !defined(__CLING__)
  typedef void (TreeWriter::*TProcessMethod)(ExRootTreeBranch *, TObjArray *); //!

  typedef std::map< ExRootTreeBranch *, std::pair< TProcessMethod, TObjArray * > > TBranchMap; //!

  TBranchMap fBranchMap; //!

  std::map< TClass *, TProcessMethod > fClassMap; //!
#endif

  ClassDef(TreeWriter, 1)
};

#endif
