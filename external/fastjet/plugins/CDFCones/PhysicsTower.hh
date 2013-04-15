#ifndef _PHYSICS_TOWER_HH_
#define _PHYSICS_TOWER_HH_

//----------------------------------------------------------------------
// This file distributed with FastJet has been obtained from
// http://www.pa.msu.edu/~huston/Les_Houches_2005/JetClu+Midpoint-StandAlone.tgz
//
// Permission to distribute it with FastJet has been granted by Joey
// Huston (see the COPYING file in the main FastJet directory for
// details).
// Changes from the original file are listed below.
//----------------------------------------------------------------------

// History of changes compared to the original PhysicsTower.hh file
// 
// 2009-01-17  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * put the code in the fastjet::cdf namespace
// 
// 2008-08-15  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * switched JetClu plugin over to proper indexed tracking of
//          jet contents rather than a (dodgy) map based on particle
//          energies 
// 
// 2006-09-24  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * replaced the private m_index variable by a public fjindex
//          one for tracking within FastJet
// 
// 2006-09-24  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * added JetClu+MidPoint to FastJet

#include "LorentzVector.hh"
#include "CalTower.hh"

#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace cdf{

class PhysicsTower
{
 public:

  LorentzVector fourVector;
  CalTower calTower;

  PhysicsTower(): fourVector(LorentzVector()), calTower(CalTower()), fjindex(-1) {}
  PhysicsTower(LorentzVector v, CalTower c): fourVector(v), calTower(c), fjindex(-1) {}
  PhysicsTower(const PhysicsTower& p): fourVector(p.fourVector), calTower(p.calTower), fjindex(p.fjindex) {}
  PhysicsTower(CalTower c):
    fourVector(LorentzVector(c.Et*cos(c.phi),c.Et*sin(c.phi),c.Et*sinh(c.eta),c.Et*cosh(c.eta))), calTower(c), fjindex(-1) {}
  PhysicsTower(LorentzVector v): fourVector(v), calTower(CalTower(v.Et(),v.eta(),v.phi())), fjindex(-1) {}
  double Et()   const {return calTower.Et;}
  double eta()  const {return calTower.eta;}
  double phi()  const {return calTower.phi;}
  int    iEta() const {return calTower.iEta;}
  int    iPhi() const {return calTower.iPhi;}
  bool isEqual(PhysicsTower p)
  {
    return fourVector.isEqual(p.fourVector) && calTower.isEqual(p.calTower);
  }
  /// addition by GPS (2008-08-15) for tracking within fastjet
  int fjindex;
};

} // namespace cdf

FASTJET_END_NAMESPACE

#endif
