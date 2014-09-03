//----------------------------------------------------------------------
// This file distributed with FastJet has been obtained from
// http://www.pa.msu.edu/~huston/Les_Houches_2005/JetClu+Midpoint-StandAlone.tgz
//
// Permission to distribute it with FastJet has been granted by Joey
// Huston (see the COPYING file in the main FastJet directory for
// details).
// Changes from the original file are listed below.
//----------------------------------------------------------------------

// History of changes compared to the original JetCluAlgorithm.cc file
// 
// 2014-08-13 Matteo Cacciari and Gavin Salam
//        * commented out towers variable in JetCluAlgorithm::buildPreClusters
//          interface to avoid compiler warning
//
// 2011-11-14  Gregory Soyez  <soyez@fastjet.fr>
//
//        * added a few parentheses suggested by the -Wparentheses gcc option
//
// 2009-01-17  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * put the code in the fastjet::cdf namespace
// 
// 2006-09-24  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * added JetClu+MidPoint to FastJet

#include "JetCluAlgorithm.hh"
#include "ClusterComparisons.hh"
#include "Centroid.hh"
#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace cdf{

void JetCluAlgorithm::makeSeedTowers(std::vector<PhysicsTower>& towers, std::vector<Cluster>& seedTowers)
{
  for(int iEta = 4; iEta < 48; iEta++){
    bool seg24 = true;
    if ((iEta >= 8 && iEta < 14) || (iEta >= 38 && iEta < 44))
      seg24 = false;
    for(int iPhi = 0; iPhi < 24; iPhi++){
      Cluster seed;
      for(std::vector<PhysicsTower>::iterator towerIter = towers.begin(); towerIter != towers.end(); towerIter++)
	if(towerIter->iEta() == iEta &&
	   ((seg24 && towerIter->iPhi() == iPhi) || (!seg24 && (towerIter->iPhi() == 2*iPhi || towerIter->iPhi() == 2*iPhi + 1))))
	  seed.addTower(*towerIter);
      if(seed.centroid.Et > _seedThreshold)
	seedTowers.push_back(seed);
    }
  }
  sort(seedTowers.begin(),seedTowers.end(),ClusterCentroidEtGreater());
}

// MC+GPS 2014-08-13, commented out the towers variable to avoid an
// unused variable warning
void JetCluAlgorithm::buildPreClusters(std::vector<Cluster>& seedTowers, std::vector<PhysicsTower>& /*towers*/,
				       std::vector<Cluster>& preClusters)
{
  std::vector<Centroid> leadingSeedTowers;
  for(std::vector<Cluster>::iterator seedTowerIter = seedTowers.begin(); seedTowerIter != seedTowers.end(); seedTowerIter++){
    bool seedTowerAddedToPreCluster = false;
    std::vector<Cluster>::iterator preClusterIter = preClusters.begin();
    std::vector<Centroid>::iterator leadingSeedTowerIter = leadingSeedTowers.begin();
    while(preClusterIter != preClusters.end() && !seedTowerAddedToPreCluster){
      double dEta = fabs(seedTowerIter->centroid.eta - leadingSeedTowerIter->eta);
      double dPhi = fabs(seedTowerIter->centroid.phi - leadingSeedTowerIter->phi);
      if(dPhi > M_PI)
	dPhi = 2*M_PI - dPhi;
      if(dEta <= _coneRadius && dPhi <= _coneRadius){
	int iEtaSeedTower = seedTowerIter->towerList.begin()->iEta();
	int iPhiSeedTower = seedTowerIter->towerList.begin()->iPhi();
	if ((iEtaSeedTower >= 8 && iEtaSeedTower < 14) || (iEtaSeedTower >= 38 && iEtaSeedTower < 44))
	  iPhiSeedTower = iPhiSeedTower/2;
 	for(std::vector<PhysicsTower>::iterator preClusterTowerIter = preClusterIter->towerList.begin();
	    preClusterTowerIter != preClusterIter->towerList.end() && !seedTowerAddedToPreCluster;
	    preClusterTowerIter++){
	  int iEtaPreClusterTower = preClusterTowerIter->iEta();
	  int iPhiPreClusterTower = preClusterTowerIter->iPhi();
	  if ((iEtaPreClusterTower >= 8 && iEtaPreClusterTower < 14) || (iEtaPreClusterTower >= 38 && iEtaPreClusterTower < 44))
	    iPhiPreClusterTower = iPhiPreClusterTower/2;
	  int dIEta = std::abs(iEtaSeedTower - iEtaPreClusterTower);
	  int dIPhi = std::abs(iPhiSeedTower - iPhiPreClusterTower);
	  if(dIPhi > 12)
	    dIPhi = 24 - dIPhi;
	  int adj = dIPhi*dIPhi + dIEta*dIEta;
	  if(adj <= _adjacencyCut){
	    for(std::vector<PhysicsTower>::iterator seedTowerTowerIter = seedTowerIter->towerList.begin();
		seedTowerTowerIter != seedTowerIter->towerList.end();
		seedTowerTowerIter++)
	      preClusterIter->addTower(*seedTowerTowerIter);
	    seedTowerAddedToPreCluster = true;
	  }
	}
      }
      preClusterIter++;
      leadingSeedTowerIter++;
    }
    if(!seedTowerAddedToPreCluster){
      Cluster newPreCluster;
      for(std::vector<PhysicsTower>::iterator seedTowerTowerIter = seedTowerIter->towerList.begin();
	  seedTowerTowerIter != seedTowerIter->towerList.end();
	  seedTowerTowerIter++)
	newPreCluster.addTower(*seedTowerTowerIter);
      preClusters.push_back(newPreCluster);
      leadingSeedTowers.push_back(Centroid(newPreCluster.centroid.Et,newPreCluster.centroid.eta,newPreCluster.centroid.phi));
    }
  }
}

void JetCluAlgorithm::findStableCones(std::vector<Cluster>& preClusters, std::vector<PhysicsTower>& towers,
				      std::vector<Cluster>& stableCones)
{
  for(std::vector<Cluster>::iterator preClusterIter = preClusters.begin(); preClusterIter != preClusters.end(); preClusterIter++){
    double startEt  = preClusterIter->centroid.Et;
    double startEta = preClusterIter->centroid.eta;
    double startPhi = preClusterIter->centroid.phi;
    int nIterations = 0;
    Cluster trialCone;
    while(nIterations++ < _maxIterations){
      trialCone.clear();
      for(std::vector<PhysicsTower>::iterator towerIter = towers.begin(); towerIter != towers.end(); towerIter++){
	double dEta = fabs(towerIter->eta() - startEta);
	double dPhi = fabs(towerIter->phi() - startPhi);
	if(dPhi > M_PI)
	  dPhi = 2*M_PI - dPhi;
	double dR = sqrt(dEta*dEta + dPhi*dPhi);
	if(dR < _coneRadius)
	  trialCone.addTower(*towerIter);
      }
      if(_iratch != 0)
	for(std::vector<PhysicsTower>::iterator preClusterTowerIter = preClusterIter->towerList.begin();
	    preClusterTowerIter != preClusterIter->towerList.end();
	    preClusterTowerIter++){
	  bool foundInTrialCone = false;
	  for(std::vector<PhysicsTower>::iterator trialConeTowerIter = trialCone.towerList.begin();
	      trialConeTowerIter != trialCone.towerList.end() && !foundInTrialCone;
	      trialConeTowerIter++)
	    if(trialConeTowerIter->isEqual(*preClusterTowerIter))
	      foundInTrialCone = true;
	  if(!foundInTrialCone)
	    trialCone.addTower(*preClusterTowerIter);
	}
      if(nIterations <= _maxIterations){
	double endEt  = trialCone.centroid.Et;
	double endEta = trialCone.centroid.eta;
	double endPhi = trialCone.centroid.phi;
	if(endEt == startEt && endEta == startEta && endPhi == startPhi)
	  nIterations = _maxIterations;
	else{
	  startEt  = endEt;
	  startEta = endEta;
	  startPhi = endPhi;
	}
      }
    }
//    bool foundIdentical = false;
//    for(std::vector<Cluster>::iterator stableConeIter = stableCones.begin();
//	stableConeIter != stableCones.end() && !foundIdentical;
//	stableConeIter++)
//      if(trialCone.centroid.isEqual(stableConeIter->centroid))
//	foundIdentical = true;
//    if(!foundIdentical)
      stableCones.push_back(trialCone);
  }
  sort(stableCones.begin(),stableCones.end(),ClusterCentroidEtGreater());
}

void JetCluAlgorithm::splitAndMerge(std::vector<Cluster>& stableCones, std::vector<Cluster>& jets)
{
  std::vector<bool> isActive;
  for(std::vector<Cluster>::iterator stableConeIter = stableCones.begin(); stableConeIter != stableCones.end(); stableConeIter++)
    isActive.push_back(bool(true));
  std::vector<bool>::iterator isActiveIter1 = isActive.begin();
  for(std::vector<Cluster>::iterator stableConeIter1 = stableCones.begin();
      stableConeIter1 != stableCones.end();
      stableConeIter1++, isActiveIter1++){
    std::vector<Cluster>::iterator stableConeIter2 = stableCones.begin();
    std::vector<bool>::iterator isActiveIter2 = isActive.begin();
    while(stableConeIter2 != stableConeIter1 && *isActiveIter1){
      if(*isActiveIter2){
	Cluster overlap;
	for(std::vector<PhysicsTower>::iterator towerIter1 = stableConeIter1->towerList.begin();
	    towerIter1 != stableConeIter1->towerList.end();
	    towerIter1++)
	  for(std::vector<PhysicsTower>::iterator towerIter2 = stableConeIter2->towerList.begin();
	      towerIter2 != stableConeIter2->towerList.end();
	      towerIter2++)
	    if(towerIter1->isEqual(*towerIter2)){
	      overlap.addTower(*towerIter1);
	      break;
	    }
	if(overlap.size()){
	  if(overlap.size() == stableConeIter2->size())
	    *isActiveIter2 = false;
	  else if(overlap.size() == stableConeIter1->size())
	    *isActiveIter1 = false;
	  else if(overlap.centroid.Et > _overlapThreshold*stableConeIter1->centroid.Et ||
		  overlap.centroid.Et > _overlapThreshold*stableConeIter2->centroid.Et){
	    for(std::vector<PhysicsTower>::iterator stableConeTowerIter2 = stableConeIter2->towerList.begin();
		stableConeTowerIter2 != stableConeIter2->towerList.end();
		stableConeTowerIter2++){
	      bool isInOverlap = false;
	      for(std::vector<PhysicsTower>::iterator overlapTowerIter = overlap.towerList.begin();
		  overlapTowerIter != overlap.towerList.end() && !isInOverlap;
		  overlapTowerIter++)
		if(stableConeTowerIter2->isEqual(*overlapTowerIter))
		  isInOverlap = true;
	      if(!isInOverlap)
		stableConeIter1->addTower(*stableConeTowerIter2);
	    }
	    *isActiveIter2 = false;
	  }
	  else{
	    Cluster removeFromStableCone1,removeFromStableCone2,oldRemoveFromStableCone1,oldRemoveFromStableCone2;
	    double etaStableCone1 = stableConeIter1->centroid.eta;
	    double phiStableCone1 = stableConeIter1->centroid.phi;
	    double etaStableCone2 = stableConeIter2->centroid.eta;
	    double phiStableCone2 = stableConeIter2->centroid.phi;
	    double dRstableCone1,dRstableCone2;
	    int iterCount = 0;
	    while(iterCount++ <= _maxIterations){
	      oldRemoveFromStableCone1.clear();
	      oldRemoveFromStableCone2.clear();
	      if(iterCount > 1){
		if(removeFromStableCone1.size()){
		  Centroid stableConeCentroid1(stableConeIter1->centroid);
		  Centroid removeCentroid1(removeFromStableCone1.centroid);
		  stableConeCentroid1.subtract(removeCentroid1);
		  etaStableCone1 = stableConeCentroid1.eta;
		  phiStableCone1 = stableConeCentroid1.phi;
		}
		else{
		  etaStableCone1 = stableConeIter1->centroid.eta;
		  phiStableCone1 = stableConeIter1->centroid.phi;
		}
		if(removeFromStableCone2.size()){
		  Centroid stableConeCentroid2(stableConeIter2->centroid);
		  Centroid removeCentroid2(removeFromStableCone2.centroid);
		  stableConeCentroid2.subtract(removeCentroid2);
		  etaStableCone2 = stableConeCentroid2.eta;
		  phiStableCone2 = stableConeCentroid2.phi;
		}
		else{
		  etaStableCone2 = stableConeIter2->centroid.eta;
		  phiStableCone2 = stableConeIter2->centroid.phi;
		}
		for(std::vector<PhysicsTower>::iterator removeTowerIter1 = removeFromStableCone1.towerList.begin();
		    removeTowerIter1 != removeFromStableCone1.towerList.end();
		    removeTowerIter1++)
		  oldRemoveFromStableCone1.addTower(*removeTowerIter1);
		for(std::vector<PhysicsTower>::iterator removeTowerIter2 = removeFromStableCone2.towerList.begin();
		    removeTowerIter2 != removeFromStableCone2.towerList.end();
		    removeTowerIter2++)
		  oldRemoveFromStableCone2.addTower(*removeTowerIter2);
	      }
	      removeFromStableCone1.clear();
	      removeFromStableCone2.clear();
	      for(std::vector<PhysicsTower>::iterator overlapTowerIter = overlap.towerList.begin();
		  overlapTowerIter != overlap.towerList.end();
		  overlapTowerIter++){
		double dEta1 = fabs(overlapTowerIter->eta() - etaStableCone1);
		double dPhi1 = fabs(overlapTowerIter->phi() - phiStableCone1);
		if(dPhi1 > M_PI)
		  dPhi1 = 2*M_PI - dPhi1;
		dRstableCone1 = dEta1*dEta1 + dPhi1*dPhi1;
		double dEta2 = fabs(overlapTowerIter->eta() - etaStableCone2);
		double dPhi2 = fabs(overlapTowerIter->phi() - phiStableCone2);
		if(dPhi2 > M_PI)
		  dPhi2 = 2*M_PI - dPhi2;
		dRstableCone2 = dEta2*dEta2 + dPhi2*dPhi2;
		if(dRstableCone1 < dRstableCone2)
		  removeFromStableCone2.addTower(*overlapTowerIter);
		else
		  removeFromStableCone1.addTower(*overlapTowerIter);
	      }
	      if(iterCount > 1 &&
		 removeFromStableCone1.size() == oldRemoveFromStableCone1.size() &&
		 removeFromStableCone2.size() == oldRemoveFromStableCone2.size() &&
		 (!removeFromStableCone1.size() || !removeFromStableCone2.size() ||
		  (removeFromStableCone1.centroid.isEqual(oldRemoveFromStableCone1.centroid) &&
		   removeFromStableCone2.centroid.isEqual(oldRemoveFromStableCone2.centroid))))
		iterCount = _maxIterations + 1;
	    }
	    for(std::vector<PhysicsTower>::iterator removeTowerIter1 = removeFromStableCone1.towerList.begin();
		removeTowerIter1 != removeFromStableCone1.towerList.end();
		removeTowerIter1++)
	      stableConeIter1->removeTower(*removeTowerIter1);
	    for(std::vector<PhysicsTower>::iterator removeTowerIter2 = removeFromStableCone2.towerList.begin();
		removeTowerIter2 != removeFromStableCone2.towerList.end();
		removeTowerIter2++)
	      stableConeIter2->removeTower(*removeTowerIter2);
	  }
	  overlap.clear();
	}
      }
      stableConeIter2++;
      isActiveIter2++;
    }
  }
  jets.clear();
  std::vector<bool>::iterator isActiveIter = isActive.begin();
  for(std::vector<Cluster>::iterator stableConeIter = stableCones.begin();
      stableConeIter != stableCones.end();
      stableConeIter++, isActiveIter++)
    if(*isActiveIter)
      jets.push_back(*stableConeIter);
  sort(jets.begin(),jets.end(),ClusterFourVectorEtGreater());
}

void JetCluAlgorithm::run(std::vector<PhysicsTower>& towers, std::vector<Cluster>& jets)
{
  std::vector<Cluster> seedTowers;
  makeSeedTowers(towers,seedTowers);
  std::vector<Cluster> preClusters;
  buildPreClusters(seedTowers,towers,preClusters);
  std::vector<Cluster> stableCones;
  findStableCones(preClusters,towers,stableCones);
  splitAndMerge(stableCones,jets);
}

}  // namespace cdf

FASTJET_END_NAMESPACE
