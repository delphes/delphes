#ifndef _JET_HH_
#define _JET_HH_

//----------------------------------------------------------------------
// This file distributed with FastJet has been obtained from SpartyJet
// v2.20.0 by Pierre-Antoine Delsart, Kurtis L. Geerlings, Joey
// Huston, Brian T. Martin and Chris Vermilion
// For details, see http://www.pa.msu.edu/~huston/SpartyJet/
//                  http://projects.hepforge.org/spartyjet/
//
// Changes from the original file are listed below.
//----------------------------------------------------------------------

// History of changes from the original Jet.hh file in SpartyJet v2.20
//  
// 2011-11-14  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * removed some harmless warnings coming with the -Wshadow gcc option
// 
// 2011-06-28  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * used stable_sort instead of sort to fix some ordering issues
// 
// 2009-01-15  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * put the code in the fastjet::atlas namespace

#include "LorentzVector.hh"
#include <list>
#include <vector>
#include <algorithm>

#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace atlas { 

class Jet : public LorentzVector {
public :
  
  typedef std::list<Jet*> constit_vect_t;
  typedef std::vector<Jet*> jet_list_t;
  
  Jet(): LorentzVector(0,0,0,0) {}
  Jet(double p1, double p2, double p3, double p0, int index_in=0): LorentzVector(p1,p2,p3,p0), m_index(index_in){}
  Jet(LorentzVector v): LorentzVector(v)  {m_index = 0;}
  Jet(Jet &j);
  Jet(Jet *j);
  
  
  /// The standard way of merging jets
  void addJet(Jet& j);
  void addJet(Jet* j);

 
  /// Access jet constituents
  int getConstituentNum(){return m_constituents.size();}
  constit_vect_t::iterator firstConstituent(){ return m_constituents.begin();};
  constit_vect_t::iterator lastConstituent(){ return m_constituents.end();};

  

  // convenience methods 
  void addConstituent(Jet* jet) {m_constituents.push_back(jet);this->add(*jet);};
  void addConstituent(constit_vect_t::iterator first, constit_vect_t::iterator last);
  void removeConstituent(Jet* jet) {m_constituents.remove(jet);this->subtract(*jet);};

  void addConstituent_notMoment(Jet* jet){m_constituents.push_back(jet);}
  //void removeConstituent(constit_vect_t::iterator first, constit_vect_t::iterator last);


  // return position in intial collection
  int index() const {return m_index;}  
  void set_index(int i){m_index= i;}


  /// Atlas compatibility code :
  LorentzVector hlv() {return *this;}


  //bool split_merged;   // from siscone/jetclu/midpoint algorithms

protected :
  int m_index;  /// position in a jet list (used for constituents positions)
  constit_vect_t m_constituents;

};


 
void find_jet_in_list(Jet* j);

// using functors is supposed to be faster... (??)
class JetSorter_Et {
public:
  bool operator()(Jet* j1, Jet* j2){
    // deal with numerical uncertainty : 
    if(fabs( j1->et() - j2->et())<0.001 ) return 0;
    else return j1->et() > j2->et();
    //return (j1->et() > j2->et());    
  }
};

class JetSorter_Pt {
public:
  bool operator()(Jet* j1, Jet* j2){
    return (j1->pt() > j2->pt());
  }
};

class JetSorter_Eta {
public:
  bool operator()(Jet* j1, Jet* j2){
    return (j1->eta() > j2->eta());
  }
};

class JetSorter_E {
public:
  bool operator()(Jet* j1, Jet* j2){
    return (j1->e() > j2->e());
  }
};



template<class T>
inline void sort_jet_list(Jet::jet_list_t &list){
  std::stable_sort(list.begin(),list.end(), T());
}
inline void sort_list_et(Jet::jet_list_t &list){
  //std::sort(list.begin(),list.end(),et_compare);
  std::stable_sort(list.begin(),list.end(), JetSorter_Et());
}
inline void sort_list_pt(Jet::jet_list_t &list){
  std::stable_sort(list.begin(),list.end(),JetSorter_Pt());
}

Jet* jet_from_overlap(Jet* j1, Jet* j2);

}  // namespace atlas

FASTJET_END_NAMESPACE
#endif


