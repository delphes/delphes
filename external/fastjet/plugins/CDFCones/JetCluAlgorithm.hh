#ifndef _JETCLU_ALGORITHM_HH_
#define _JETCLU_ALGORITHM_HH_

//----------------------------------------------------------------------
// This file distributed with FastJet has been obtained from
// http://www.pa.msu.edu/~huston/Les_Houches_2005/JetClu+Midpoint-StandAlone.tgz
//
// Permission to distribute it with FastJet has been granted by Joey
// Huston (see the COPYING file in the main FastJet directory for
// details).
// Changes from the original file are listed below.
//----------------------------------------------------------------------

// History of changes compared to the original JetCluAlgorithm.hh file
// 
// 2009-01-17  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * put the code in the fastjet::cdf namespace
// 
// 2006-09-24  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * added JetClu+MidPoint to FastJet

#include "PhysicsTower.hh"
#include "Cluster.hh"
#include <vector>

#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace cdf{

class JetCluAlgorithm
{
 private:
  double _seedThreshold;
  double _coneRadius;
  int    _adjacencyCut;
  int    _maxIterations;
  int    _iratch;
  double _overlapThreshold;

 public:
  JetCluAlgorithm():
    _seedThreshold(1),
    _coneRadius(0.7),
    _adjacencyCut(2),
    _maxIterations(100),
    _iratch(1),
    _overlapThreshold(0.75)
  {}
  JetCluAlgorithm(double st, double cr, int ac, int mi, int ir, double ot):
    _seedThreshold(st),
    _coneRadius(cr),
    _adjacencyCut(ac),
    _maxIterations(mi),
    _iratch(ir),
    _overlapThreshold(ot)
  {}
  void makeSeedTowers(std::vector<PhysicsTower>& towers, std::vector<Cluster>& seedTowers);
  void buildPreClusters(std::vector<Cluster>& seedTowers, std::vector<PhysicsTower>& towers, std::vector<Cluster>& preClusters);
  void findStableCones(std::vector<Cluster>& preClusters, std::vector<PhysicsTower>& towers, std::vector<Cluster>& stableCones);
  void splitAndMerge(std::vector<Cluster>& stableCones, std::vector<Cluster>& jets);
  void run(std::vector<PhysicsTower>& towers, std::vector<Cluster>& jets);
};

}  // namespace cdf

FASTJET_END_NAMESPACE

#endif
