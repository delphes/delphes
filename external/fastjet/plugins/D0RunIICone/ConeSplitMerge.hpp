#ifndef  D0RunIIconeJets_CONESPLITMERGE
#define  D0RunIIconeJets_CONESPLITMERGE
// ---------------------------------------------------------------------------
// ConeSplitMerge.hpp
//
// Created: 28-JUL-2000 Francois Touze
//
// Purpose: Implements the pT ordered jet split/merge algorithm for the 
//   Improved Legacy Cone Algorithm split/merge algo.
//
// Modified:
//   31-JUL-2000  Laurent Duflot
//     + introduce support for additional informations (ConeJetInfo)
//    1-AUG-2000  Laurent Duflot
//     + jet merge criteria was wrong, now calculate shared_ET and compare to 
//       neighbour jet pT.
//     + split was incomplete: shared items not really removed from jets.
//    4-Aug-2000  Laurent Duflot
//     + use item methods y() and phi() rather than p4vec() and then P2y and 
//       P2phi
//    7-Aug-2000  Laurent Duflot 
//     + force the list to be organized by decreasing ET and, for identical ET,
//       by decreasing seedET. Identical protojets can be found by eg nearby
//       seeds. The seedET ordering is such that for identical jets, the one
//       with the highest seedET is kept, which is what we want for efficiency
//       calculations.
//     + avoid unecessary copies of lists by using reference when possible
//    9-Aug-2000  Laurent Duflot
//     + save initial jet ET before split/merge
//    1-May-2007 Lars Sonnenschein
//    extracted from D0 software framework and modified to remove subsequent dependencies
//
//
// This file is distributed with FastJet under the terms of the GNU
// General Public License (v2). Permission to do so has been granted
// by Lars Sonnenschein and the D0 collaboration (see COPYING for
// details)
//
// History of changes in FastJet compared tothe original version of
// ConeSplitMerge.hpp
//
// 2011-12-13  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * added license information
//
// 2009-01-17  Gregory Soyez  <soyez@fastjet.fr>
//
//        * put the code in the fastjet::d0 namespace
//
// 2007-12-14  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * replaced make_pair by std::make_pair
//
// ---------------------------------------------------------------------------


#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include "ProtoJet.hpp"

//using namespace D0RunIIconeJets_CONEJETINFO;

#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace d0{

//
// this class is used to order ProtoJets by decreasing ET and seed ET
template <class Item>
class ProtoJet_ET_seedET_order
{
public:
  bool operator()(const ProtoJet<Item> & first, const ProtoJet<Item> & second) const
  {
    if ( first.pT() > second.pT() ) return true;
    else
      if ( first.pT() < second.pT() ) return false;
      else return ( first.info().seedET() > second.info().seedET() );
  }
};


template <class Item>
class ConeSplitMerge {

public :

  typedef std::multimap<ProtoJet<Item>,float,ProtoJet_ET_seedET_order<Item> > PJMMAP;
  
  ConeSplitMerge();
  ConeSplitMerge(const std::vector<ProtoJet<Item> >& jvector);
  ConeSplitMerge(const std::list<ProtoJet<Item> >& jlist);
  ~ConeSplitMerge() {;}
  void split_merge(std::vector<ProtoJet<Item> >& ecv,float s, float pT_min_leading_protojet, float pT_min_second_protojet, int MergeMax, float pT_min_noMergeMax);

private :
  PJMMAP _members;
};
///////////////////////////////////////////////////////////////////////////////
template<class Item>
ConeSplitMerge<Item>::ConeSplitMerge():_members() {;}

template<class Item>
ConeSplitMerge<Item>::ConeSplitMerge(const std::vector<ProtoJet<Item> >& jvector) 
{
  // sort proto_jets in Et descending order
  typename std::vector<ProtoJet<Item> >::const_iterator jt;
  for(jt = jvector.begin(); jt != jvector.end(); ++jt) 
  {
    // this is supposed to be a stable cone, declare so
    ProtoJet<Item> jet(*jt);
    jet.NowStable();
    _members.insert(std::make_pair(jet,jet.pT()));
  }
}

template<class Item>
ConeSplitMerge<Item>::ConeSplitMerge(const std::list<ProtoJet<Item> >& jlist) 
{
  //_max_nb_items =-1;
  // sort proto_jets in Et descending order
  typename std::list<ProtoJet<Item> >::const_iterator jt;
  for(jt = jlist.begin(); jt != jlist.end(); ++jt) 
  {
    // this is supposed to be a stable cone, declare so
    ProtoJet<Item> jet(*jt);
    jet.NowStable();
    _members.insert(std::make_pair(jet,jet.pT()));
  }
}

template<class Item>
void ConeSplitMerge<Item>::split_merge(std::vector<ProtoJet<Item> >& jcv,
				       float shared_ET_fraction,
				       float pT_min_leading_protojet, float pT_min_second_protojet,
				       int MergeMax, float pT_min_noMergeMax) 
{
  while(!_members.empty()) 
  {
    /*
    {
      std::cout << std::endl;
      std::cout << " ---------------  list of protojets ------------------ " <<std::endl;
      for ( PJMMAP::iterator it = _members.begin();
	    it != _members.end(); ++it)
      {
	std::cout << " pT y phi " << (*it).pT() << " " << (*it).y() << " " << (*it).phi() << " " << (*it).info().seedET() <<  " " << (*it).info().nbMerge() << " " << (*it).info().nbSplit() << std::endl;
      }
      std::cout << " ----------------------------------------------------- " <<std::endl;
    }
    */


    // select highest Et jet
    typename PJMMAP::iterator itmax= _members.begin();
    ProtoJet<Item> imax((*itmax).first);
    const std::list<const Item*>& ilist(imax.LItems());

    _members.erase(itmax);
 
    // does jet share items?
    bool share= false;
    float shared_ET = 0.;
    typename PJMMAP::iterator jtmax;
    typename PJMMAP::iterator jt;
    for(jt = _members.begin(); jt != _members.end(); ++jt) 
    {
      const std::list<const Item*>& jlist((*jt).first.LItems());
      typename std::list<const Item*>::const_iterator tk;
      int count;
      for(tk = ilist.begin(), count = 0; tk != ilist.end(); 
	  ++tk, ++count) 
      {
	typename std::list<const Item*>::const_iterator where= 
	  find(jlist.begin(),jlist.end(),*tk);   
	if(where != jlist.end()) 
	{
	  share= true;
	  shared_ET += (*tk)->pT();
	}
      }
      if(share) 
      {
	jtmax = jt;
	break;
      }
    }
    
    if(!share) {
      // add jet to the final jet list
      jcv.push_back(imax);
      //std::cout << " final jet " << imax.pT() << " "<< imax.info().nbMerge() << " " << imax.info().nbSplit() << std::endl; 
    }
    else 
    {
      // find highest Et neighbor
      ProtoJet<Item> jmax((*jtmax).first);

      // drop neighbor
      _members.erase(jtmax);


      //std::cout << " split or merge ? " << imax.pT() << " " << shared_ET << " " << jmax.pT() << " " << s << " " << (jmax.pT())*s << std::endl;

      // merge
      if( shared_ET > (jmax.pT())*shared_ET_fraction 
	  && (imax.pT() > pT_min_leading_protojet || jmax.pT() > pT_min_second_protojet)
	  && (imax.info().nbMerge() < MergeMax || imax.pT() > pT_min_noMergeMax))
	{
	  // add neighbor's items to imax
	  const std::list<const Item*>& jlist(jmax.LItems());
	  typename std::list<const Item*>::const_iterator tk;
	  typename std::list<const Item*>::const_iterator iend= ilist.end();
	  bool same = true; // is jmax just the same as imax ? 
	  for(tk = jlist.begin(); tk != jlist.end(); ++tk) 
	    {
	      typename std::list<const Item*>::const_iterator where= 
		find(ilist.begin(),iend,*tk);   
	      if(where == iend) 
		{
		  imax.addItem(*tk);
		  same = false;
		}
	    }
	  if ( ! same ) 
	    {
	      // recalculate
	      //float old_pT = imax.pT();
	      
	      imax.updateJet();
	      imax.merged();
	      //std::cout << " jet merge :: " << old_pT << " " << jmax.pT() << " " << imax.pT() << " "<< imax.info().nbMerge() << " " << imax.info().nbSplit() << std::endl; 
	    }
	}
      
      //split and assign removed shared cells from lowest pT protojet
      else if(shared_ET > (jmax.pT())*shared_ET_fraction)
      {
	// here we need to pull the lists, because there are items to remove                                                                           
	std::list<const Item*> ilist(imax.LItems());
	std::list<const Item*> jlist(jmax.LItems());

        typename std::list<const Item*>::iterator tk;
        for(tk = jlist.begin(); tk != jlist.end(); )
	  {
	    typename std::list<const Item*>::iterator where=
	      find(ilist.begin(),ilist.end(),*tk);
	    if(where != ilist.end()) {
	      tk = jlist.erase(tk);
	    }
	    else ++tk;
	  }
	
        jmax.erase();

        for ( typename std::list<const Item*>::const_iterator it = jlist.begin();
              it != jlist.end(); ++it) jmax.addItem(*it);

        // recalculated jet quantities 
        jmax.updateJet();
        jmax.splitted();
        //std::cout << " jet split 1 :: " << jmax.pT() << " "<< jmax.info().nbMerge() << " " << jmax.info().nbSplit() << std::endl;                         
        _members.insert(std::make_pair(jmax,jmax.pT()));
      }

      // split and assign shared cells to nearest protojet
      else 
      {
	// here we need to pull the lists, because there are items to remove
	std::list<const Item*> ilist(imax.LItems());
	std::list<const Item*> jlist(jmax.LItems());
	
	typename std::list<const Item*>::iterator tk;
	for(tk = jlist.begin(); tk != jlist.end(); ) 
	{
	  typename std::list<const Item*>::iterator where= 
	    find(ilist.begin(),ilist.end(),*tk);
	  if(where != ilist.end()) {
	    float yk   = (*tk)->y();
	    float phik = (*tk)->phi();
	    float di= RD2(imax.y(),imax.phi(),yk,phik);
	    float dj= RD2(jmax.y(),jmax.phi(),yk,phik);
	    if(dj > di) 
	    {
	      tk = jlist.erase(tk);
	      //std::cout << " shared item " << (*tk)->pT() << " removed from neighbour jet " << std::endl;
	    }
	    else
	    {
	      ilist.erase(where);
	      ++tk;
	      //std::cout << " shared item " << (*tk)->pT() << " removed from leading jet " << std::endl;
	    }
	  }
	  else ++tk;
	}
	// recalculate jets imax and jmax
	
	// first erase all items
	imax.erase();
	// put items that are left
	for ( typename std::list<const Item*>::const_iterator it = ilist.begin();
	      it != ilist.end(); ++it) imax.addItem(*it);
	// recalculated jet quantities
	imax.updateJet();
	imax.splitted();
	//std::cout << " jet split 2 :: " << imax.pT() << " "<< imax.info().nbMerge() << " " << imax.info().nbSplit() << std::endl; 


	// first erase all items
	jmax.erase();
	// put items that are left
	for ( typename std::list<const Item*>::const_iterator it = jlist.begin();
	      it != jlist.end(); ++it) jmax.addItem(*it);
	// recalculated jet quantities
	jmax.updateJet();
	jmax.splitted();
	//std::cout << " jet split " << jmax.pT() << " "<< jmax.info().nbMerge() << " " << jmax.info().nbSplit() << std::endl; 

	_members.insert(std::make_pair(jmax,jmax.pT()));
      }
      _members.insert(std::make_pair(imax,imax.pT()));
    }
  } // while
}
///////////////////////////////////////////////////////////////////////////////

}  // namespace d0

FASTJET_END_NAMESPACE

#endif
