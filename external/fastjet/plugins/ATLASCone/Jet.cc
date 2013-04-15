//----------------------------------------------------------------------
// This file distributed with FastJet has been obtained from SpartyJet
// v2.20.0 by Pierre-Antoine Delsart, Kurtis L. Geerlings, Joey
// Huston, Brian T. Martin and Chris Vermilion
// For details, see http://www.pa.msu.edu/~huston/SpartyJet/
//                  http://projects.hepforge.org/spartyjet/
//
// Changes from the original file are listed below.
//----------------------------------------------------------------------

// History of changes from the original Jet.cc file in SpartyJet v2.20
//  
// 2009-01-15  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * put the code in the fastjet::atlas namespace

#include "Jet.hh"
#include <iostream>

#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace atlas { 

Jet::Jet(Jet &jet) : LorentzVector(0,0,0,0){
  add(jet);
  m_index = jet.index();
  m_constituents = jet.m_constituents;
  //  m_area = jet.area();
  //  m_area_error = jet.area_error();
}
Jet::Jet(Jet* j){
  add(*j);
  m_index = j->index();
  m_constituents = j->m_constituents;
//   m_area = j->area();
//   m_area_error = j->area_error();
}

void Jet::addJet(Jet& j){
  add(j);
  m_constituents.insert(m_constituents.end(), j.firstConstituent(), j.lastConstituent() );
}

void Jet::addJet(Jet* j){
  add(*j) ;
  m_constituents.insert(m_constituents.end(), j->firstConstituent(), j->lastConstituent() );
}


void Jet::addConstituent(constit_vect_t::iterator first, constit_vect_t::iterator last){
  m_constituents.insert(m_constituents.end(), first, last); 
  for(; first!=last;++first) this->add( **first);
}



Jet* jet_from_overlap(Jet* j1, Jet* j2){
  Jet *j = new Jet();
  Jet::constit_vect_t::iterator it1 = j1->firstConstituent();
  Jet::constit_vect_t::iterator it1E = j1->lastConstituent();
  for(;it1!= it1E; ++it1){
    Jet::constit_vect_t::iterator it2 = j2->firstConstituent();
    Jet::constit_vect_t::iterator it2E = j2->lastConstituent();
    for(;it2!= it2E; ++it2){
      if( *it1 == *it2) j->addConstituent(*it1);
    }
    
  }
  return j;
}
}  // namespace atlas

FASTJET_END_NAMESPACE

