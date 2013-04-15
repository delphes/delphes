//----------------------------------------------------------------------
// This file distributed with FastJet has been obtained from
// http://www.pa.msu.edu/~huston/Les_Houches_2005/JetClu+Midpoint-StandAlone.tgz
//
// Permission to distribute it with FastJet has been granted by Joey
// Huston (see the COPYING file in the main FastJet directory for
// details).
// Changes from the original file are listed below.
//----------------------------------------------------------------------

// History of changes compared to the original MidPointAlgorithm.cc file
// 
// 2009-01-17  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * put the code in the fastjet::cdf namespace
// 
// 2008-01-15  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * put in the volatile fix for midpoint with optimization;
//          checked regression test, and that speed comes out OK
// 
// 
// 2007-10-16  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * added protection to midpoint for cases where no stable cones
//          are found
// 
// 2007-03-10  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * added support for the pttilde scale choice in the CDF midpoint code
//
// 2007-02-21  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * added option of choosing the scale used in the split-merge
//          procedure (pt [default], Et or mt)
// 
// 2006-09-24  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * removed the local_sort method
// 
// 2006-09-24  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * added JetClu+MidPoint to FastJet

#include "MidPointAlgorithm.hh"
#include "ClusterComparisons.hh"
#include <algorithm>
#include <iostream>
#include <cmath>

#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace cdf{

void MidPointAlgorithm::findStableConesFromSeeds(std::vector<PhysicsTower>& towers, std::vector<Cluster>& stableCones)
{
  bool reduceConeSize = true;
  for(std::vector<PhysicsTower>::iterator towerIter = towers.begin(); towerIter != towers.end(); towerIter++)
    if(towerIter->fourVector.pt() > _seedThreshold)
      iterateCone(towerIter->fourVector.y(),towerIter->fourVector.phi(),0,towers,stableCones,reduceConeSize);
}

void MidPointAlgorithm::findStableConesFromMidPoints(std::vector<PhysicsTower>& towers, std::vector<Cluster>& stableCones)
{
  // distanceOK[i-1][j] = Is distance between stableCones i and j (i>j) less than 2*_coneRadius?
  std::vector< std::vector<bool> > distanceOK;
  distanceOK.resize(stableCones.size() - 1);
  for(int nCluster1 = 1; nCluster1 < stableCones.size(); nCluster1++){
    distanceOK[nCluster1 - 1].resize(nCluster1);
    double cluster1Rapidity = stableCones[nCluster1].fourVector.y();
    double cluster1Phi      = stableCones[nCluster1].fourVector.phi();
    for(int nCluster2 = 0; nCluster2 < nCluster1; nCluster2++){
      double cluster2Rapidity = stableCones[nCluster2].fourVector.y();
      double cluster2Phi      = stableCones[nCluster2].fourVector.phi();
      double dRapidity = fabs(cluster1Rapidity - cluster2Rapidity);
      double dPhi      = fabs(cluster1Phi      - cluster2Phi);
      if(dPhi > M_PI)
        dPhi = 2*M_PI - dPhi;
      double dR = sqrt(dRapidity*dRapidity + dPhi*dPhi);
      distanceOK[nCluster1 - 1][nCluster2] = dR < 2*_coneRadius;
    }
  }

  // Find all pairs (triplets, ...) of stableCones which are less than 2*_coneRadius apart from each other.
  std::vector< std::vector<int> > pairs(0);
  std::vector<int> testPair(0);
  int maxClustersInPair = _maxPairSize;
  if(!maxClustersInPair)
    maxClustersInPair = stableCones.size();
  addClustersToPairs(testPair,pairs,distanceOK,maxClustersInPair);

  // Loop over all combinations. Calculate MidPoint. Make midPointClusters.
  bool reduceConeSize = false;
  for(int iPair = 0; iPair < pairs.size(); iPair++){
    // Calculate rapidity, phi and pT of MidPoint.
    LorentzVector midPoint(0,0,0,0);
    for(int iPairMember = 0; iPairMember < pairs[iPair].size(); iPairMember++)
      midPoint.add(stableCones[pairs[iPair][iPairMember]].fourVector);
    iterateCone(midPoint.y(),midPoint.phi(),midPoint.pt(),towers,stableCones,reduceConeSize);
  }

  //sort(stableCones.begin(),stableCones.end(),ClusterPtGreater());
  local_sort(stableCones);  // GPS mod. to allow split-merge with various scales
}


void MidPointAlgorithm::iterateCone(volatile double startRapidity, volatile double startPhi, volatile double startPt,
				    std::vector<PhysicsTower>& towers, std::vector<Cluster>& stableCones, bool reduceConeSize)
{
  int nIterations = 0;
  bool keepJet = true;
  Cluster trialCone;
  double iterationConeRadius = _coneRadius;
  if(reduceConeSize)
    iterationConeRadius *= sqrt(_coneAreaFraction);
  while(nIterations++ < _maxIterations + 1 && keepJet){
    trialCone.clear();
    // Find particles which should go in the cone.
    if(nIterations == _maxIterations + 1)
      iterationConeRadius = _coneRadius;
    for(std::vector<PhysicsTower>::iterator towerIter = towers.begin(); towerIter != towers.end(); towerIter++){
      double dRapidity = fabs(towerIter->fourVector.y()   - startRapidity);
      double dPhi      = fabs(towerIter->fourVector.phi() - startPhi);
      if(dPhi > M_PI)
	dPhi = 2*M_PI - dPhi;
      double dR = sqrt(dRapidity*dRapidity + dPhi*dPhi);
      if(dR < iterationConeRadius)
	trialCone.addTower(*towerIter);
    }
    if(!trialCone.size())   // Empty cone?
      keepJet = false;
    else{
      if(nIterations <= _maxIterations){
	volatile double endRapidity = trialCone.fourVector.y();
	volatile double endPhi      = trialCone.fourVector.phi();
	volatile double endPt       = trialCone.fourVector.pt();
	// Do we have a stable cone?
	if(endRapidity == startRapidity && endPhi == startPhi && endPt == startPt){
	  // If cone size is reduced, then do one more iteration.
	  nIterations = _maxIterations;
	  if(!reduceConeSize)
	    nIterations++;
	}
	else{
	  // Another iteration.
	  startRapidity = endRapidity;
	  startPhi      = endPhi;
	  startPt       = endPt;
	}
      }
    }
  }

  if(keepJet){
    // We have a stable cone.
    bool identical = false;
    for(std::vector<Cluster>::iterator stableConeIter = stableCones.begin(); stableConeIter != stableCones.end(); stableConeIter++)
      if(trialCone.fourVector.isEqual(stableConeIter->fourVector))
	identical = true;
    if(!identical)
      stableCones.push_back(trialCone);
  }
}

void MidPointAlgorithm::addClustersToPairs(std::vector<int>& testPair, std::vector< std::vector<int> >& pairs,
					   std::vector< std::vector<bool> >& distanceOK, int maxClustersInPair)
{
  // Recursively adds clusters to pairs, triplets, ... whose mid-points are then calculated.

  // Find StableCone number to start with (either 0 at the beginning or last element of testPair + 1).
  int nextClusterStart = 0;
  if(testPair.size())
    nextClusterStart = testPair.back() + 1;
  for(int nextCluster = nextClusterStart; nextCluster <= distanceOK.size(); nextCluster++){
    // Is new SeedCone less than 2*_coneRadius apart from all clusters in testPair?
    bool addCluster = true;
    for(int iCluster = 0; iCluster < testPair.size() && addCluster; iCluster++)
      if(!distanceOK[nextCluster - 1][testPair[iCluster]])
	addCluster = false;
    if(addCluster){
      // Add it to the testPair.
      testPair.push_back(nextCluster);
      // If testPair is a pair, add it to pairs.
      if(testPair.size() > 1)
	pairs.push_back(testPair);
      // If not bigger than allowed, find more clusters within 2*_coneRadius.
      if(testPair.size() < maxClustersInPair)
	addClustersToPairs(testPair,pairs,distanceOK,maxClustersInPair);
      // All combinations containing testPair found. Remove last element.
      testPair.pop_back();
    }
  }
}

void MidPointAlgorithm::splitAndMerge(std::vector<Cluster>& stableCones, std::vector<Cluster>& jets)
{
  bool mergingNotFinished = true;
  while(mergingNotFinished){
    // Sort the stable cones (highest pt first).
    //sort(stableCones.begin(),stableCones.end(),ClusterPtGreater());
    local_sort(stableCones);

    // Start with the highest pt cone.
    std::vector<Cluster>::iterator stableConeIter1 = stableCones.begin();
    if(stableConeIter1 == stableCones.end())   // Stable cone list empty?
      mergingNotFinished = false;
    else{
      bool coneNotModified = true;
      // Determine whether highest pt cone has an overlap with other stable cones.
      std::vector<Cluster>::iterator stableConeIter2 = stableConeIter1;
      stableConeIter2++;   // 2nd highest pt cone.
      while(coneNotModified && stableConeIter2 != stableCones.end()){
	// Calculate overlap of the two cones.
	Cluster overlap;
	for(std::vector<PhysicsTower>::iterator towerIter1 = stableConeIter1->towerList.begin();
	    towerIter1 != stableConeIter1->towerList.end();
	    towerIter1++){
	  bool isInCone2 = false;
	  for(std::vector<PhysicsTower>::iterator towerIter2 = stableConeIter2->towerList.begin();
	      towerIter2 != stableConeIter2->towerList.end();
	      towerIter2++)
	    if(towerIter1->isEqual(*towerIter2))
	      isInCone2 = true;
	  if(isInCone2)
	    overlap.addTower(*towerIter1);
	}
	if(overlap.size()){   // non-empty overlap
	  coneNotModified = false;
          // GPS mod to allow various scale choices in split merge --------
          double overlap_scale, jet2_scale;
          switch(_smScale) {
          case SM_pt:
            overlap_scale = overlap.fourVector.pt();
            jet2_scale    = stableConeIter2->fourVector.pt();
            break;
          case SM_Et:
            overlap_scale = overlap.fourVector.Et();
            jet2_scale    = stableConeIter2->fourVector.Et();
            break;
          case SM_mt:
            overlap_scale = overlap.fourVector.mt();
            jet2_scale    = stableConeIter2->fourVector.mt();
            break;
          case SM_pttilde:
            overlap_scale = overlap.pt_tilde;
            jet2_scale    = stableConeIter2->pt_tilde;
            break;
          default:
            std::cerr << "Unrecognized value for _smScale: " 
                      << _smScale << std::endl;
            exit(-1);
          }
	  if(overlap_scale >= _overlapThreshold*jet2_scale){
          // end of GPS modification ---------------------------
          //if(overlap.fourVector.pt() >= _overlapThreshold*stableConeIter2->fourVector.pt()){
	    // Merge the two cones.
	    for(std::vector<PhysicsTower>::iterator towerIter2 = stableConeIter2->towerList.begin();
		towerIter2 != stableConeIter2->towerList.end();
		towerIter2++){
	      bool isInOverlap = false;
	      for(std::vector<PhysicsTower>::iterator overlapTowerIter = overlap.towerList.begin();
		  overlapTowerIter != overlap.towerList.end();
		  overlapTowerIter++)
		if(towerIter2->isEqual(*overlapTowerIter))
		  isInOverlap = true;
	      if(!isInOverlap)
		stableConeIter1->addTower(*towerIter2);
	    }
	    // Remove the second cone.
	    stableCones.erase(stableConeIter2);
	  }
	  else{
	    // Separate the two cones.
	    // Which particle goes where?
	    std::vector<PhysicsTower> removeFromCone1,removeFromCone2;
	    for(std::vector<PhysicsTower>::iterator towerIter = overlap.towerList.begin();
		towerIter != overlap.towerList.end();
		towerIter++){
	      double towerRapidity = towerIter->fourVector.y();
	      double towerPhi      = towerIter->fourVector.phi();
	      // Calculate distance from cone 1.
	      double dRapidity1 = fabs(towerRapidity - stableConeIter1->fourVector.y());
	      double dPhi1      = fabs(towerPhi      - stableConeIter1->fourVector.phi());
	      if(dPhi1 > M_PI)
		dPhi1 = 2*M_PI - dPhi1;
	      double dRJet1 = sqrt(dRapidity1*dRapidity1 + dPhi1*dPhi1);
	      // Calculate distance from cone 2.
	      double dRapidity2 = fabs(towerRapidity - stableConeIter2->fourVector.y());
	      double dPhi2      = fabs(towerPhi      - stableConeIter2->fourVector.phi());
	      if(dPhi2 > M_PI)
		dPhi2 = 2*M_PI - dPhi2;
	      double dRJet2 = sqrt(dRapidity2*dRapidity2 + dPhi2*dPhi2);
	      if(dRJet1 < dRJet2)
		// Particle is closer to cone 1. To be removed from cone 2.
		removeFromCone2.push_back(*towerIter);
	      else
		// Particle is closer to cone 2. To be removed from cone 1.
		removeFromCone1.push_back(*towerIter);
	    }
	    // Remove particles in the overlap region from the cones to which they have the larger distance.
	    for(std::vector<PhysicsTower>::iterator towerIter = removeFromCone1.begin();
		towerIter != removeFromCone1.end();
		towerIter++)
	      stableConeIter1->removeTower(*towerIter);
	    for(std::vector<PhysicsTower>::iterator towerIter = removeFromCone2.begin();
		towerIter != removeFromCone2.end();
		towerIter++)
	      stableConeIter2->removeTower(*towerIter);
	  }
	}
	stableConeIter2++;
      }
      if(coneNotModified){
	// Cone 1 has no overlap with any of the other cones and can become a jet.
	jets.push_back(*stableConeIter1);
	stableCones.erase(stableConeIter1);
      }
    }
  }

  //sort(jets.begin(),jets.end(),ClusterPtGreater());
  local_sort(jets); // GPS mod. to allow split-merge with various scales
}



void MidPointAlgorithm::local_sort(std::vector<Cluster>& clusters) {
  switch(_smScale) {
  case SM_pt:
    sort(clusters.begin(),clusters.end(),ClusterPtGreater());
    break;
  case SM_Et:
    sort(clusters.begin(),clusters.end(),ClusterFourVectorEtGreater());
    break;
  case SM_mt:
    sort(clusters.begin(),clusters.end(),ClusterMtGreater());
    break;
  case SM_pttilde:
    sort(clusters.begin(),clusters.end(),ClusterPtTildeGreater());
    break;
  default:
    std::cerr << "Unrecognized value for _smScale: " << _smScale << std::endl;
    exit(-1);
  }
}



void MidPointAlgorithm::run(std::vector<PhysicsTower>& towers, std::vector<Cluster>& jets)
{
  std::vector<Cluster> stableCones;
  findStableConesFromSeeds(towers,stableCones);
  // GPS addition to prevent crashes if no stable cones
  // are found (e.g. all particles below seed threshold)
  if (stableCones.size() > 0) {
    findStableConesFromMidPoints(towers,stableCones);
    splitAndMerge(stableCones,jets);
  }
}

}  // namespace cdf

FASTJET_END_NAMESPACE
