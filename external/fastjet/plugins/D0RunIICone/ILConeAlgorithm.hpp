#ifndef  D0RunIIconeJets_ILCONEALGORITHM
#define  D0RunIIconeJets_ILCONEALGORITHM
// ---------------------------------------------------------------------------
// ILConeAlgorithm.hpp
//
// Created: 28-JUL-2000 Francois Touze (+ Laurent Duflot)
//
// Purpose: Implements the Improved Legacy Cone Algorithm
//
// Modified:
//   31-JUL-2000  Laurent Duflot
//     + introduce support for additional informations (ConeJetInfo)
//    1-AUG-2000  Laurent Duflot
//     + seedET for midpoint jet is -999999., consistent with seedET ordering
//       in ConeSplitMerge: final jets with seedET=-999999. will be midpoint 
//       jets which were actually different from stable cones.
//    4-Aug-2000  Laurent Duflot
//     + remove unecessary copy in TemporaryJet::is_stable()
//   11-Aug-2000  Laurent Duflot
//     + remove using namespace std
//     + add threshold on item. The input list to makeClusters() IS modified
//   20-June-2002 John Krane
//     + remove bug in midpoint calculation based on pT weight
//     + started to add new midpoint calculation based on 4-vectors,
//       but not enough info held by this class
//   24-June-2002 John Krane
//     + modify is_stable() to not iterate if desired 
//       (to expand search cone w/out moving it)
//     + added search cone logic for initial jets but not midpoint jets as per
//       agreement with CDF
//   19-July-2002 John Krane
//     + _SEARCH_CONE size factor now provided by calreco/CalClusterReco.cpp 
//     + (default = 1.0 = like original ILCone behavior)
//   10-Oct-2002 John Krane
//     + Check min Pt of cone with full size first, then iterate with search cone
//   07-Dec-2002 John Krane
//     + speed up the midpoint phi-wrap solution
//   01-May-2007 Lars Sonnenschein
//   extracted from D0 software framework and modified to remove subsequent dependencies
//
//
// This file is distributed with FastJet under the terms of the GNU
// General Public License (v2). Permission to do so has been granted
// by Lars Sonnenschein and the D0 collaboration (see COPYING for
// details)
//
// History of changes in FastJet compared tothe original version of
// ProtoJet.hpp
//
// 2012-06-12  Gregory Soyez  <soyez@fastjet.fr>
//        * Replaced addItem(...) by this->addItem(...) to allow
//          compilation with gcc 4.7 which no longer performs
//          unqualified template lookups. See
//          e.g. http://gcc.gnu.org/gcc-4.7/porting_to.html
//
// 2011-12-13  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * added license information
//
// 2011-11-14  Gregory Soyez  <soyez@fastjet.fr>
//
//        * changed the name of a few parameters to avoid a gcc
//          -Wshadow warning
//
// 2009-01-17  Gregory Soyez  <soyez@fastjet.fr>
//
//        * put the code in the fastjet::d0 namespace
//
// 2007-12-14  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * moved the 'std::vector<ProtoJet<Item> > ilcv' structure
//          containing the info about the final jets from a local
//          variable to a class variable (for integration in the
//          FastJet plugin core)
//
// ---------------------------------------------------------------------------

#include <vector>
#include <list>
#include <utility>  // defines pair<,>
#include <map>
#include <algorithm>
#include <iostream>


//#include "energycluster/EnergyClusterReco.hpp"
//#include "energycluster/ProtoJet.hpp"
#include "ProtoJet.hpp"
//#include "energycluster/ConeSplitMerge.hpp"
#include "ConeSplitMerge.hpp"
//#include "energycluster/ConeJetInfoChunk.hpp"

#include "inline_maths.h"

///////////////////////////////////////////////////////////////////////////////
#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace d0{

using namespace inline_maths;

/*
 NB: Some attempt at optimizing the code has been made by ordering the object
     along rapidity but this did not improve the speed. There are traces of 
     these atemps in the code, that will be cleaned up in the future.
 */

// at most one of those !
// order the input list and use lower_bound() and upper_bound() to restrict the
// on item to those that could possibly be in the cone.
//#define ILCA_USE_ORDERED_LIST

// idem but use an intermediate multimap in hope that lower_bound() and 
// upper_bound() are faster in this case.
//#define ILCA_USE_MMAP


#ifdef ILCA_USE_ORDERED_LIST
// this class is used to order the item list along rapidity
template <class Item>
class rapidity_order
{
public:
  bool operator()(const Item* first, const Item* second)
  {
    return (first->y() < second->y());
  }
  bool operator()(float const & first, const Item* second)
  {
    return (first  < second->y());
  }
  bool operator()(const Item* first, float const& second)
  {
    return (first->y() < second);
  }
};
#endif


//template <class Item,class ItemAddress,class IChunk>
template <class Item>
class ILConeAlgorithm 
{

public:

  // default constructor (default parameters are crazy: you should not use that
  // constructor !)
  ILConeAlgorithm():
    _CONE_RADIUS(0.),
    _MIN_JET_ET(99999.),
    _ET_MIN_RATIO(0.5),
    _FAR_RATIO(0.5),
    _SPLIT_RATIO(0.5),
    _DUPLICATE_DR(0.005),
    _DUPLICATE_DPT(0.01),
    _SEARCH_CONE(0.5),
    _PT_MIN_LEADING_PROTOJET(0), 
    _PT_MIN_SECOND_PROTOJET(0), 
    _MERGE_MAX(10000), 
    _PT_MIN_noMERGE_MAX(0)
  {;}

  // full constructor
  ILConeAlgorithm(float cone_radius, float min_jet_Et, float split_ratio,
		  float far_ratio=0.5, float Et_min_ratio=0.5,
		  bool kill_duplicate=true, float duplicate_dR=0.005, 
		  float duplicate_dPT=0.01, float search_factor=1.0, 
		  float pT_min_leading_protojet=0., float pT_min_second_protojet=0.,
		  int merge_max=10000, float pT_min_nomerge=0.) :
    // cone radius
    _CONE_RADIUS(cone_radius), 
    // minimum jet ET
    _MIN_JET_ET(min_jet_Et), 
    // stable cones must have ET > _ET_MIN_RATIO*_MIN_JET_ET at any iteration
    _ET_MIN_RATIO(Et_min_ratio),
    // precluster at least _FAR_RATIO*_CONE_RADIUS away from stable cones
    _FAR_RATIO(far_ratio), 
    // split or merge criterium           
    _SPLIT_RATIO(split_ratio),
    _DUPLICATE_DR(duplicate_dR),
    _DUPLICATE_DPT(duplicate_dPT),
    _SEARCH_CONE(cone_radius/search_factor),
    // kill stable cone if within _DUPLICATE_DR and delta(pT)<_DUPLICATE_DPT
    // of another stable cone.
    _KILL_DUPLICATE(kill_duplicate),
    _PT_MIN_LEADING_PROTOJET(pT_min_leading_protojet),
    _PT_MIN_SECOND_PROTOJET(pT_min_second_protojet),
    _MERGE_MAX(merge_max),
    _PT_MIN_noMERGE_MAX(pT_min_nomerge)
    {;}

  //destructor
  ~ILConeAlgorithm() {;}

  // make jet clusters using improved legacy cone algorithm
  //void makeClusters(const EnergyClusterReco* r,
  void makeClusters(
		    // the input list modified (ordered)
		    std::list<Item> &jets,
		    std::list<const Item*>& itemlist, 
		    //float zvertex,   
		    ////std::list<const Item*>& itemlist);   
		    //const EnergyClusterCollection<ItemAddress>& preclu,
		    //IChunk* chunkptr, ConeJetInfoChunk* infochunkptr,
		    const float Item_ET_Threshold);

  // this will hold the final jets + contents
  std::vector<ProtoJet<Item> > ilcv;

private:

  float _CONE_RADIUS;
  float _MIN_JET_ET;
  float _ET_MIN_RATIO;
  float _FAR_RATIO;
  float _SPLIT_RATIO;
  float _DUPLICATE_DR;
  float _DUPLICATE_DPT;
  float _SEARCH_CONE;

  bool _KILL_DUPLICATE;

  float _PT_MIN_LEADING_PROTOJET;
  float _PT_MIN_SECOND_PROTOJET;
  int  _MERGE_MAX; 
  float _PT_MIN_noMERGE_MAX;

  // private class 
  // This is a ProtoJet with additional methods dist(), midpoint() and 
  // is_stable()
  class TemporaryJet : public ProtoJet<Item> 
  {
    
  public :
    
    TemporaryJet(float seedET) : ProtoJet<Item>(seedET) {;}

    TemporaryJet(float seedET,float y_in,float phi_in) : 
      ProtoJet<Item>(seedET,y_in,phi_in) {;}
    
    ~TemporaryJet() {;}
    
    float dist(TemporaryJet& jet) const 
    {
      return RDelta(this->_y,this->_phi,jet.y(),jet.phi()); 
    }
    
    void midpoint(const TemporaryJet& jet,float & y_out, float & phi_out) const 
    {
      // Midpoint should probably be computed w/4-vectors but don't 
      // have that info.  Preserving Pt-weighted calculation - JPK
      float pTsum = this->_pT + jet.pT();
      y_out = (this->_y*this->_pT + jet.y()*jet.pT())/pTsum;

      phi_out = (this->_phi*this->_pT + jet.phi()*jet.pT())/pTsum;
      // careful with phi-wrap area: convert from [0,2pi] to [-pi,pi]
      //ls: original D0 code, as of 23/Mar/2007
      //if ( abs(phi-this->_phi)>2.0 ) { // assumes cones R=1.14 or smaller, merge within 2R only  
      //ls: abs bug fixed 26/Mar/2007 
      if ( fabs(phi_out-this->_phi)>2.0 ) { // assumes cones R=1.14 or smaller, merge within 2R only  
        phi_out = fmod( this->_phi+PI, TWOPI);
	if (phi_out < 0.0) phi_out += TWOPI;
	phi_out -= PI;

	float temp=fmod( jet.phi()+PI, TWOPI);
	if (temp < 0.0) temp += TWOPI;
	temp -= PI;

	phi_out = (phi_out*this->_pT + temp*jet.pT()) /pTsum;
      }

      if ( phi_out < 0. ) phi_out += TWOPI;
    }
    

////////////////////////////////////////
#ifdef ILCA_USE_MMAP
    bool is_stable(const std::multimap<float,const Item*>& items, 
		   float radius, float min_ET, int max_iterations=50) 
#else
    bool is_stable(const std::list<const Item*>& itemlist, float radius, 
		 float min_ET, int max_iterations=50) 
#endif
    // Note: max_iterations = 0 will just recompute the jet using the specified cone
    {
      float radius2 = radius*radius;
      float Rcut= 1.E-06;
      
      
      // ?? if(_Increase_Delta_R) Rcut= 1.E-04;
      bool stable= true;
      int trial= 0;
      float Yst;
      float PHIst;
      do {  
	trial++;
	//std::cout << "   trial " << trial << " " << _y << " " << _phi << std::endl; 
	Yst  = this->_y;
	PHIst= this->_phi;    
	//cout << "is_stable beginning do loop: this->_pT=" << this->_pT << " this->_y=" << this->_y << " this->_phi=" << this->_phi << endl;
	this->erase();
	
	this->setJet(Yst,PHIst,0.0);
	
	
#ifdef ILCA_USE_ORDERED_LIST      
	std::list<const Item*>::const_iterator lower = 
	  lower_bound(itemlist.begin(),itemlist.end(),Yst-radius,
		      rapidity_order<Item>());
	std::list<const Item*>::const_iterator upper = 
	  upper_bound(itemlist.begin(),itemlist.end(),Yst+radius,
		      rapidity_order<Item>());
	for(std::list<const Item*>::const_iterator tk = lower; tk != upper; ++tk)      {
	  //std::cout << " is_stable: item y=" << (*tk)->y() << " phi=" << (*tk)->phi() << " RD2=" << RD2((*tk)->y(),(*tk)->phi(),Yst,PHIst) << " " << Yst-radius << " " << Yst+radius << endl;
	  if(RD2((*tk)->y(),(*tk)->phi(),Yst,PHIst) <= radius2) 
	    {
	      this->addItem(*tk);
	    }
	}       
#else
#ifdef ILCA_USE_MMAP      
	// need to loop only on the subset with   Yst-R < y < Yst+R
	for ( std::multimap<float,const Item*>::const_iterator 
		tk = items.lower_bound(Yst-radius);
	      tk != items.upper_bound(Yst+radius); ++tk )
	  {
	    //std::cout << "     item " << (*tk)->y() << " " << (*tk)->phi() << " " << RD2((*tk)->y(),(*tk)->phi(),Yst,PHIst) << " " << Yst-radius << " " << Yst+radius << endl;
	    if(RD2(((*tk).second)->y(),((*tk).second)->phi(),Yst,PHIst) <= radius2) 
	      {
		this->addItem((*tk).second);
	      }
	  }
	
#else   

	//cout << " is_stable: itemlist.size()=" << itemlist.size() << endl;
	for(typename std::list<const Item*>::const_iterator tk = itemlist.begin(); tk != itemlist.end(); ++tk) 
	  {
	    //std::cout << "    is_stable: item (*tk)->y()=" << (*tk)->y() << " (*tk)->phi=" << (*tk)->phi() << " RD2=" << RD2((*tk)->y(),(*tk)->phi(),Yst,PHIst) << " Yst-rad=" << Yst-radius << " Yst+rad=" << Yst+radius << endl;
	    if(RD2((*tk)->y(),(*tk)->phi(),Yst,PHIst) <= radius2) 
	       {
		 //cout << "add item to *tk" << endl;
		this->addItem(*tk);
	      }
	  }
#endif
#endif      
      
	//std::cout << "is_stable, before update: jet this->_y=" << this->_y << " _phi=" << this->_phi << " _pT=" << this->_pT << " min_ET=" << min_ET << std::endl; 
	this->updateJet();
	//std::cout << "is_stable, after update: jet this->_y=" << this->_y << " _phi=" << this->_phi << " _pT=" << this->_pT << " min_ET=" << min_ET << std::endl; 
	
	if(this->_pT < min_ET ) 
	  {
	    stable= false;
	    break;
	  } 
	//cout << "is_stable end while loop: this->_pT=" << this->_pT << " this->_y=" << this->_y << " this->_phi=" << this->_phi << endl;
      } while(RD2(this->_y,this->_phi,Yst,PHIst) >= Rcut && trial <= max_iterations);
      //std::cout << "   trial stable " << trial << " " << stable << std::endl; 
      return stable;
    }
    
  private :
    
  };
};
///////////////////////////////////////////////////////////////////////////////
//template <class Item,class ItemAddress,class IChunk>
//void ILConeAlgorithm <Item,ItemAddress,IChunk>::
template <class Item>
void ILConeAlgorithm <Item>::
//makeClusters(const EnergyClusterReco* r,
makeClusters(
             std::list<Item> &jets,
	     std::list<const Item*>& ilist, 
	     //float Zvertex, 
	     ////std::list<const Item*>& ilist) 
	     //const EnergyClusterCollection<ItemAddress>& preclu,
	     //IChunk* chunkptr, ConeJetInfoChunk* infochunkptr,
	     const float Item_ET_Threshold) 
{
  // remove items below threshold
  for ( typename std::list<const Item*>::iterator it = ilist.begin(); 

        it != ilist.end(); )
  {
    if ( (*it)->pT() < Item_ET_Threshold ) 
    {
      it = ilist.erase(it);
    }
      else ++it;
  }

  // create an energy cluster collection for jets 
  //EnergyClusterCollection<ItemAddress>* ptrcol;
  //Item* ptrcol;
  //r->createClusterCollection(chunkptr,ptrcol);
  ////std::vector<const EnergyCluster<ItemAddress>*> ecv;
  std::vector<const Item*> ecv;
  for ( typename std::list<const Item*>::iterator it = ilist.begin(); 
        it != ilist.end(); it++) {
    ecv.push_back(*it);
  }


  //preclu.getClusters(ecv);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% need to fill here vector ecv

  //std::cout << " number of seed clusters: " << ecv.size() << std::endl;

  // skip precluster close to jets
  float far_def = _FAR_RATIO*_CONE_RADIUS * _FAR_RATIO*_CONE_RADIUS;

  // skip if jet Et is below some value
  float ratio= _MIN_JET_ET*_ET_MIN_RATIO;


#ifdef ILCA_USE_ORDERED_LIST
  // sort the list in rapidity order
  ilist.sort(rapidity_order<Item>());
#else
#ifdef ILCA_USE_MMAP  
  // create a y ordered list of items 
  std::multimap<float,const Item*> items;
  std::list<const Item*>::const_iterator it;
  for(it = ilist.begin(); it != ilist.end(); ++it) 
  {
    pair<float,const Item*> p((*it)->y(),*it);
    items.insert(p);
  }
#endif
#endif

  std::vector<ProtoJet<Item> > mcoll;
  std::vector<TemporaryJet> scoll; 


  // find stable jets around seeds 
  //typename std::vector<const EnergyCluster<ItemAddress>* >::iterator jclu;
  typename std::vector<const Item*>::iterator jclu;
  for(jclu = ecv.begin(); jclu != ecv.end(); ++jclu) 
  {
    //const EnergyCluster<ItemAddress>* ptr= *jclu;
    const Item* ptr= *jclu;
    float p[4];
    ptr->p4vec(p);
    float Yst  = P2y(p);
    float PHIst= P2phi(p);

    // don't keep preclusters close to jet
    bool is_far= true;
    // ?? if(_Kill_Far_Clusters) {
    for(unsigned int i = 0; i < scoll.size(); ++i) 
    {
      float y  = scoll[i].y();
      float phi= scoll[i].phi();
      if(RD2(Yst,PHIst,y,phi) < far_def) 
      {
	is_far= false;
	break;
      }
    }
    // ?? }

    if(is_far) 
    {
      TemporaryJet jet(ptr->pT(),Yst,PHIst);
      //cout << "temporary jet: pt=" << ptr->pT() << " y=" << Yst << " phi=" << PHIst << endl;

      // Search cones are smaller, so contain less jet Et 
      // Don't throw out too many little jets during search phase!
      // Strategy: check Et first with full cone, then search with low-GeV min_et thresh
#ifdef ILCA_USE_MMAP
      if(jet.is_stable(items,_CONE_RADIUS,ratio,0) && jet.is_stable(items,_SEARCH_CONE,3.0)) 
#else
      if(jet.is_stable(ilist,_CONE_RADIUS,ratio,0) && jet.is_stable(ilist,_SEARCH_CONE,3.0)) 
#endif
      {

	//cout << "temporary jet is stable" << endl;

// jpk  Resize the found jets 
#ifdef ILCA_USE_MMAP
        jet.is_stable(items,_CONE_RADIUS,ratio,0) ;
#else
        jet.is_stable(ilist,_CONE_RADIUS,ratio,0) ;
#endif
	//cout << "found jet has been resized if any" << endl;

	if ( _KILL_DUPLICATE ) 
	{
	  // check if we are not finding the same jet again
	  float distmax = 999.; int imax = -1;
	  for(unsigned int i = 0; i < scoll.size(); ++i) 
	  {
	    float dist = jet.dist(scoll[i]);
	    if ( dist < distmax )
	    {
	      distmax = dist;
	      imax = i;
	    }
	  }
	  if ( distmax > _DUPLICATE_DR ||
	       fabs((jet.pT()-scoll[imax].pT())/scoll[imax].pT())>_DUPLICATE_DPT )
	  {
	    scoll.push_back(jet);
	    mcoll.push_back(jet);
	    //std::cout << " stable cone " << jet.y() << " " << jet.phi() << " " << jet.pT() << std::endl;
	  }
	  /*
	    else
	    {
	    std::cout << " jet too close to a found jet " << distmax << " " << 
	    fabs((jet.pT()-scoll[imax].pT())/scoll[imax].pT()) << std::endl;
	    }
	  */
	}
	else
	{
	  scoll.push_back(jet);
	  mcoll.push_back(jet);
	  //std::cout << " stable cone " << jet.y() << " " << jet.phi() << " " << jet.pT() << std::endl;
	}

      }
    }
  }

  // find stable jets around midpoints
  for(unsigned int i = 0; i < scoll.size(); ++i) 
  {
    for(unsigned int k = i+1; k < scoll.size(); ++k) 
    {
      float djet= scoll[i].dist(scoll[k]);
      if(djet > _CONE_RADIUS && djet < 2.*_CONE_RADIUS) 
      {
	float y_mid,phi_mid;
	scoll[i].midpoint(scoll[k],y_mid,phi_mid);
	TemporaryJet jet(-999999.,y_mid,phi_mid);
	//std::cout << " midpoint: " << scoll[i].pT() << " " << scoll[i].info().seedET() << " " << scoll[k].pT() << " " << scoll[k].info().seedET() << " " << y_mid << " " << phi_mid << std::endl; 

// midpoint jets are full size
#ifdef ILCA_USE_MMAP
      if(jet.is_stable(items,_CONE_RADIUS,ratio)) 
#else
      if(jet.is_stable(ilist,_CONE_RADIUS,ratio)) 
#endif
	{
	  mcoll.push_back(jet);
	  //std::cout << " stable midpoint cone " << jet.y() << " " << jet.phi() << " " << jet.pT() << std::endl;
	}
      }
    }
  }


  // do a pT ordered splitting/merging
  ConeSplitMerge<Item> pjets(mcoll);
  ilcv.clear();
  pjets.split_merge(ilcv,_SPLIT_RATIO, _PT_MIN_LEADING_PROTOJET, _PT_MIN_SECOND_PROTOJET, _MERGE_MAX, _PT_MIN_noMERGE_MAX);


  for(unsigned int i = 0; i < ilcv.size(); ++i) 
  {
    if ( ilcv[i].pT() > _MIN_JET_ET )
    {
      //EnergyCluster<ItemAddress>* ptrclu;
      Item ptrclu;
      //r->createCluster(ptrcol,ptrclu);
      
      std::list<const Item*> tlist=ilcv[i].LItems();
      typename std::list<const Item*>::iterator tk;
      for(tk = tlist.begin(); tk != tlist.end(); ++tk) 
      {
	//ItemAddress addr= (*tk)->address();
	float pk[4];
	(*tk)->p4vec(pk);
        //std::cout << (*tk)->index <<" " <<  (*tk) << std::endl;
	//float emE= (*tk)->emE();
	//r->addClusterItem(ptrclu,addr,pk,emE);
	//ptrclu->Add(*tk);
	ptrclu.Add(**tk);
      }
      // print out
      //ptrclu->print(cout);
      
      //infochunkptr->addInfo(ilcv[i].info());
      jets.push_back(ptrclu);
    }
  }
}

}  // namespace d0

FASTJET_END_NAMESPACE

#endif


