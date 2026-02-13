/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TreeWriter_h
#define TreeWriter_h

/** \class TreeWriter
 *
 *  Fills ROOT tree branches.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <map>

class TClass;
class TRefArray;

class Candidate;
class ExRootTreeBranch;

class TreeWriter : public DelphesModule
{
public:
  TreeWriter();
  ~TreeWriter();

  void Init();
  void Process();
  void Finish();

private:
  void FillParticles(const Candidate &candidate, TRefArray *array);

  void ProcessParticles(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessVertices(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessTracks(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessTowers(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessParticleFlowCandidates(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessPhotons(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessElectrons(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessMuons(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessCscCluster(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessTauJets(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessJets(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessMissingET(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessScalarHT(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessRho(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessWeight(ExRootTreeBranch *branch, const std::vector<Candidate> &array);
  void ProcessHectorHit(ExRootTreeBranch *branch, const std::vector<Candidate> &array);

#if !defined(__CINT__) && !defined(__CLING__)
  typedef void (TreeWriter::*TProcessMethod)(ExRootTreeBranch *, const std::vector<Candidate> &); //!

  typedef std::map<ExRootTreeBranch *, std::pair<TProcessMethod, const std::vector<Candidate> &> > TBranchMap; //!

  TBranchMap fBranchMap; //!

  std::map<TClass *, TProcessMethod> fClassMap; //!
#endif

  ClassDef(TreeWriter, 2)
};

#endif
