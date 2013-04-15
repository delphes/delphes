//----------------------------------------------------------------------
// This file distributed with FastJet has been obtained from SpartyJet
// v2.20.0 by Pierre-Antoine Delsart, Kurtis L. Geerlings, Joey
// Huston, Brian T. Martin and Chris Vermilion
// For details, see http://www.pa.msu.edu/~huston/SpartyJet/
//                  http://projects.hepforge.org/spartyjet/
//
// Changes from the original file are listed below.
//----------------------------------------------------------------------

//*******************************************************************************
// Filename : JetSplitMergeTool.cxx 
// Author   : Ambreesh Gupta
// Created  : Nov, 2001
//
// File taken from SpartyJet v2.20.0
// Modifications:
//   removed the string name in the ctor
//   removed the Message m_log
//   replaced px() -> px, ... in the LorentzVector calls
//   cleaned the comments
//*******************************************************************************

// History of changes from the original JetSplitMergeTool.cc file in
// SpartyJet v2.20
// 
// 2009-01-15  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * put the code in the fastjet::atlas namespace
// 
// 2009-02-14  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * imported into FastJet
//        * removed the string name in the ctor
//        * removed the Message m_log
//        * replaced px() -> px, ... in the LorentzVector calls
//        * cleaned the comments

#include <iostream>

#include "JetSplitMergeTool.hh"
#include "Jet.hh"
#include "JetDistances.hh"
#include "CommonUtils.hh"

#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace atlas { 

JetSplitMergeTool::JetSplitMergeTool()
  :  m_f( 0.5 )
{}

JetSplitMergeTool::~JetSplitMergeTool()
{}

/////////////////////////////////////////////////////////////////////////////////
//Execution                                                                     /
/////////////////////////////////////////////////////////////////////////////////
int JetSplitMergeTool::execute( jetcollection_t* theJets )
{
  m_ctr = 0;
  m_dctr = 0;

  ////////////////////////////////////////////////////
  // From the input, collection create a list of Jet//
  ////////////////////////////////////////////////////
  m_preJet.clear();
  m_jet.clear();

  jetcollection_t::iterator itrB = theJets->begin();
  jetcollection_t::iterator itrE = theJets->end(); 

  double etot =0.;
  for (;itrB!=itrE;itrB++) {
    Jet* j = new Jet(); j->addJet(*itrB);
    m_ctr +=1;
    m_preJet.push_back(j);    

    etot += j->e();    
  }

  /////////////////////
  // Split Merge Jets//
  /////////////////////
  this->split_merge();
 
  /////////////////////////////////////////////
  // Empty and re-fill input jetcollection_t //
  /////////////////////////////////////////////
  clear_list(*theJets);
  jetcollection_t::iterator it = m_jet.begin();
  jetcollection_t::iterator itE = m_jet.end();
  for(; it!=itE; ++it){    
    theJets->push_back(*it);
  }

  return 1;
}

///////////////////////////////////////////////////////////////////////////////
// Reconstruction algorithm specific methods                                  /
///////////////////////////////////////////////////////////////////////////////

void JetSplitMergeTool::split_merge()
{
  if ( m_preJet.size() >= 2 ) {
    do {
      sort_list_et(m_preJet);
      
      jetcollection_t::iterator itr;
      jetcollection_t::iterator first = m_preJet.begin();
      jetcollection_t::iterator last  = m_preJet.end();
      
      itr=first;
      ++itr;
      bool overlap = false;
  
      for (;itr != last;++itr) {      
	double etaF = (*first)->eta();
	double phiF = (*first)->phi();
	double etaS = (*itr)->eta();
	double phiS = (*itr)->phi();
	
	Jet* oJet = jet_from_overlap( (*first),*itr);   
	m_ctr +=1; 

	Jet::constit_vect_t::iterator itro  = oJet->firstConstituent();
	Jet::constit_vect_t::iterator itroE = oJet->lastConstituent();
	
	if ( oJet->getConstituentNum() != 0 ) {
	  overlap = true;
	  
	  // fraction
	  double f = sqrt(pow(oJet->px,2)+pow(oJet->py,2))/
	    sqrt(pow((*itr)->px,2)+pow((*itr)->py,2));
	  
	  // merge
	  if ( f > m_f) {
	    // we need to remove constituents !
	    Jet *j = (*first);
	    for ( ;itro != itroE; ++itro ) j->removeConstituent(*itro);
	    (*first)->addJet(*itr);
	    //m_preJet.remove(*itr);
	    delete *itr;
	    m_preJet.erase(itr);
	    m_dctr +=1;
	  }	
	  
	  // split	
	  if ( f <= m_f) {	  
	    for ( ;itro != itroE; ++itro ) {	      
	      // Distance of first jet from ProtoJet
	      double deta1 = etaF - (*itro)->eta();
	      double dphi1 = fabs(JetDistances::deltaPhi(phiF,(*itro)->phi()));
	      double dist1 = pow( deta1 , 2 ) + pow( dphi1 , 2 );
	      
	      // Distance of second jet from ProtoJet
	      double deta2 = etaS - (*itro)->eta();
	      double dphi2 = fabs(JetDistances::deltaPhi(phiS,(*itro)->phi()));
	      double dist2 = pow( deta2 , 2 ) + pow( dphi2 , 2 );
	      
	      // Remove protojet from farther Jet	    	      
	      if ( dist1 > dist2 ) (*first)->removeConstituent(*itro);
	      if ( dist1 <= dist2 ) (*itr)->removeConstituent(*itro);	
	    }
	  }
	  // Delete overlap jet     
	  delete oJet;     
	  m_dctr +=1;
	  break; 
	}  
	else {
	  // Delete overlap jet     
	  delete oJet;     
	  m_dctr +=1;
	}    
      }
      
      if ( overlap == false ) {
	m_jet.push_back(*first);
	//m_preJet.remove(*first);      
	m_preJet.erase(first);      
      }
      
    } while ( m_preJet.size() != 0 );    
  }
  else if ( m_preJet.size() == 1) {
    m_jet.push_back( *(m_preJet.begin()) );
  }

}

//////////////////////////////////////////////////////////////////////

// "True" eta and phi ASSUMING the 4-vector is filled as
// ex -> e * sin(theta) * cos(phi)
// ey -> e * sin(theta) * sin(phi)
// ez -> e * cos(theta)
// e  -> e 
// Jet phi range is (-pi,pi].


double JetSplitMergeTool::etaTrue(Jet::constit_vect_t::iterator pj)
{
  double s = ((*pj)->e() > 0) ? +1.0 : -1.0;
  double px = (*pj)->px;
  double py = (*pj)->py;
  double pz = (*pj)->pz;
  double theta = acos(pz*s/sqrt(px*px+py*py+pz*pz));
  return -log(tan(theta/2.0));
}

double JetSplitMergeTool::phiTrue(Jet::constit_vect_t::iterator pj)
{
  double s = ((*pj)->e() > 0) ? +1.0 : -1.0;
  double px = (*pj)->px;
  double py = (*pj)->py;
  return atan2(py*s,px*s);
}

}  // namespace atlas

FASTJET_END_NAMESPACE
