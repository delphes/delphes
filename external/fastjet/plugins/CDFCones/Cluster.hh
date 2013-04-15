#ifndef _CLUSTER_HH_
#define _CLUSTER_HH_

//----------------------------------------------------------------------
// This file distributed with FastJet has been obtained from
// http://www.pa.msu.edu/~huston/Les_Houches_2005/JetClu+Midpoint-StandAlone.tgz
//
// Permission to distribute it with FastJet has been granted by Joey
// Huston (see the COPYING file in the main FastJet directory for
// details).
// Changes from the original file are listed below.
//----------------------------------------------------------------------

// History of changes compared to the original Cluster.hh file
// 
// 2009-01-17  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * put the code in the fastjet::cdf namespace
// 
// 2007-03-10  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * added support for the pttilde scale choice in the CDF midpoint code
// 
// 2006-09-24  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * added JetClu+MidPoint to FastJet

#include "PhysicsTower.hh"
#include "LorentzVector.hh"
#include "Centroid.hh"
#include <vector>

#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace cdf{

class Cluster
{
 public:
  std::vector<PhysicsTower> towerList;
  LorentzVector fourVector;
  Centroid centroid;
  // addition by G.P.Salam; pt_tilde = sum |p_{ti}|. Maintaining this
  // seems to add about 1% (3%) to overall timings for midpoint
  // (jetclu) but it is useful because it makes it easy to look at
  // other scales in the split-merge procedure
  double pt_tilde; 

  Cluster()
  {
    clear();
  }
  void clear()
  {
    towerList.clear();
    fourVector = LorentzVector();
    centroid = Centroid();
    pt_tilde = 0.0;
  }
  void addTower(PhysicsTower p)
  {
    towerList.push_back(p);
    fourVector.add(p.fourVector);
    centroid.add(Centroid(p.Et(),p.eta(),p.phi()));
    pt_tilde += p.fourVector.pt();
  }
  void removeTower(PhysicsTower p)
  {
    for(std::vector<PhysicsTower>::iterator towerIter = towerList.begin(); towerIter != towerList.end(); towerIter++)
      if(towerIter->isEqual(p)){
	fourVector.subtract(towerIter->fourVector);
	centroid.subtract(Centroid(towerIter->Et(),towerIter->eta(),towerIter->phi()));
        pt_tilde -= towerIter->fourVector.pt();
	towerList.erase(towerIter);
	break;
      }
  }
  int size(){return towerList.size();}
};

}  // namespace cdf

FASTJET_END_NAMESPACE

#endif
