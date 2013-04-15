//////////////////////////////////////////////////////////////
//  File: ConeClusterAlgo.hpp                                 
//
//  Author: G. Le Meur & F. Touze
//
//  Created: 15-JUNE-1998        
//
//  Purpose: make jet clusters using fixed cone like algorithm
//           implemented in RUNI.
//
//  Modified: 
//  28-OCT-1998 to use KinemUtil (S. Protopopescu)
//   8-JAN-1999: Laurent Duflot
//     . correct bugs in getItemsInCone and updateEtaPhiEt for jets
//       overlapping the phi=0 line
//     . change abs(float) to fabs(float)
//   1-NOV-1999: Laurent Duflot
//     . correct bug in makeCluster: when the temporary jet was emptied the eta
//       and phi were not set again. The main effect was a nearly zero 
//       efficiency for jets at phi=pi (as seen by Volker Buescher)
//   25-JAN-2000: Francois Touze
//     . change in updateEtaPhiEt : the method E() which returns energy doesn't
//       exist in MCparticle classe,... so use the 4-momentum components
//     . declare const the EnergyClusterCollection of seeds in makeClusters
//   01-FEB-2000: Laurent Duflot
//     . add a missing break statement in the removal of shared items. Caused
//       an infinite loop on some events.
//     . correct typo in variable name. Change a variable name that was in 
//       French.
//     . leave some debug printout (commented)
//   15-Sep-2009 Lars Sonnenschein
//   extracted from D0 software framework and modified to remove subsequent dependencies
//
//
// This file is distributed with FastJet under the terms of the GNU
// General Public License (v2). Permission to do so has been granted
// by Lars Sonnenschein and the D0 collaboration (see COPYING for
// details)
//
// History of changes in FastJet compared to the original version of
// ConeClusterAlgo.hpp
//
// 2011-12-13  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * added license information
//
// 2011-10-06  Gregory Soyez  <soyez@fastjet.fr>
//
//        * put the code in the fastjet::d0runi namespace
//
//////////////////////////////////////////////////////////////

//#ifndef CONECLUSTERALGO_H
//#define CONECLUSTERALGO_H

#ifndef  D0RunIconeJets_CONECLUSTERALGO_H
#define  D0RunIconeJets_CONECLUSTERALGO_H


//#include "EnergyClusterReco.hpp"
#include <vector>
#include <list>
#include <utility>
//#include "kinem_util/AnglesUtil.hpp"

#include <algorithm>
#include <iostream>

#include "inline_maths.h"
#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace d0runi{

using namespace std;

//some utility functions
inline float R2(float eta1, float phi1, float eta2, float phi2) {
  return (eta1-eta2)*(eta1-eta2)+(phi1-phi2)*(phi1-phi2); }

inline float R2_bis(float eta1, float phi1, float eta2, float phi2) {
  //float dphi = kinem::delta_phi(phi1,phi2);
  float dphi = inline_maths::delta_phi(phi1,phi2);
  return (eta1-eta2)*(eta1-eta2)+dphi*dphi; }

inline float DELTA_r(float eta1,float eta2,float phi1,float phi2) {
  //float dphi = kinem::delta_phi(phi1,phi2);
  float dphi = inline_maths::delta_phi(phi1,phi2);
  return sqrt((eta1-eta2)*(eta1-eta2)+dphi*dphi);
}

inline float E2eta(float* p) { 
   float small= 1.E-05;
   float E[3];
   if(p[3] < 0.0) {
    E[0]= -p[0];
    E[1]= -p[1];
    E[2]= -p[2];
   }
   else {
    E[0]= p[0];
    E[1]= p[1];
    E[2]= p[2];
   }
   float pperp= sqrt(E[0]*E[0]+E[1]*E[1])+small;
   float ptotal= sqrt(E[0]*E[0]+E[1]*E[1]+E[2]*E[2])+small;
   //float theta= atan2(pperp,E[2]);
 
   float eta= 0.0;
   if(E[2] > 0.0) eta= log((ptotal+E[2])/pperp);
   else eta= log(pperp/(ptotal-E[2]));
   return eta;
}

inline float E2phi(float* p) { 
   float small= 1.E-05;
   float E[3];
   if(p[3] < 0.0) {
    E[0]= -p[0];
    E[1]= -p[1];
    E[2]= -p[2];
   }
   else {
    E[0]= p[0];
    E[1]= p[1];
    E[2]= p[2];
   }
   float phi= atan2(E[1],E[0]+small);
   //if(phi < 0.0) phi+=kinem::TWOPI;
   if (phi < 0.0) phi += inline_maths::TWOPI;
   return phi;
}

//template < class CalItem,class CalItemAddress,class CalIClusterChunk >
template < class CalItem >
class ConeClusterAlgo {
  //
  // Purpose: make calorimeter clusters using a cone algorithm from 
  // preclusters created previously by the class ConePreClusterAlgo.
  // Items must have addresses and 4-momenta.
  // The algorithm is implemented with a template function makeClusters.
  //  
  public :


//default constructor
ConeClusterAlgo() {} 

//constructor for cone jet algorithm
ConeClusterAlgo( float CONErad,float JETmne,float SPLifr):
  _CONErad(fabs(CONErad)), 
  _JETmne(JETmne), 
  _SPLifr(SPLifr),
  _TWOrad(0.),
  _D0_Angle(false),
  _Increase_Delta_R(true),
  _Kill_Far_Clusters(true),
  _Jet_Et_Min_On_Iter(true),
  _Far_Ratio(0.5),
  _Eitem_Negdrop(-1.0),
  _Et_Min_Ratio(0.5),
  _Thresh_Diff_Et(0.01)
  {}

//changing default thresholds & parameters
// (declared by PARAMETER in RUNI code)
ConeClusterAlgo( float CONErad,float JETmne,float SPLifr,float TWOrad, 
                 float Tresh_Diff_Et,bool D0_Angle,bool Increase_Delta_R,
                 bool Kill_Far_Clusters,bool Jet_Et_Min_On_Iter,
                 float Far_Ratio,float Eitem_Negdrop,float Et_Min_Ratio ):
  _CONErad(fabs(CONErad)), 
  _JETmne(JETmne),
  _SPLifr(SPLifr),
  _TWOrad(TWOrad),
  _D0_Angle(D0_Angle),
  _Increase_Delta_R(Increase_Delta_R),
  _Kill_Far_Clusters(Kill_Far_Clusters),
  _Jet_Et_Min_On_Iter(Jet_Et_Min_On_Iter),
  _Far_Ratio(Far_Ratio),
  _Eitem_Negdrop(Eitem_Negdrop),
  _Et_Min_Ratio(Et_Min_Ratio),
  _Thresh_Diff_Et(Tresh_Diff_Et)
  {}

//destructor
~ConeClusterAlgo() {}
  
//to make jet clusters using cone algorithm
void makeClusters(//const EnergyClusterReco* r,
		  std::list<CalItem> &jets,
		  list<const CalItem*> &itemlist, float Zvertex 
		  //, const EnergyClusterCollection<CalItemAddress> &preclu,
		  //CalIClusterChunk* chunkptr
		  //) const;
		  );

//print parameters of the algorithm
void print(ostream &os)const;

  //vector< TemporaryJet > TempColl;  


  private :
  
  float _CONErad;
  float _JETmne;
  float _SPLifr;

  float _TWOrad;
  bool _D0_Angle;
  bool _Increase_Delta_R;
  bool _Kill_Far_Clusters;
  bool _Jet_Et_Min_On_Iter;
  float _Far_Ratio;
  float _Eitem_Negdrop;
  float _Et_Min_Ratio;
  float _Thresh_Diff_Et;

  class TemporaryJet {

  public:



    TemporaryJet() {}

    TemporaryJet(float eta,float phi) { 
      _Eta=eta; 
      _Phi=phi;
    }

    ~TemporaryJet() {}

    void addItem(const CalItem* tw) {
      _LItems.push_back(tw);
    }

    void setEtaPhiEt(float eta,float phi,float pT) {
      _Eta= eta;
      _Phi= phi;
      _Et = pT;
    }

    void erase() {
      _LItems.erase(_LItems.begin(),_LItems.end());
      _Eta= 0.0;
      _Phi= 0.0;
      _Et = 0.0;  
    }

    bool share_jets(TemporaryJet &NewJet,float SharedFr,float SPLifr) {
      //
      // combined
      //
      if(SharedFr >= SPLifr) {
	typename list<const CalItem*>::iterator it;
	typename list<const CalItem*>::iterator end_of_old=_LItems.end();
	for(it=NewJet._LItems.begin(); it!=NewJet._LItems.end(); it++) {
	  typename list<const CalItem*>::iterator 
	    where = find(_LItems.begin(),end_of_old,*it);
	  // if the item is not shared, add to this jet
	  if(where == end_of_old) {
            _LItems.push_back(*it);
          }
	}
	NewJet.erase();
	return false;
      } else {
	//
	// split
	//
	typename list<const CalItem*>::iterator it;
	for(it=NewJet._LItems.begin(); it!=NewJet._LItems.end(); ) {
	  typename list<const CalItem*>::iterator 
	    where = find(_LItems.begin(),_LItems.end(),*it);
	  if(where != _LItems.end()) {
	    //float EtaItem=(*it)->eta();
	    //float PhiItem=(*it)->phi();
	    // stay closer to the RUNI conventions for negative E cells
	    float pz[4];
	    (*it)->p4vec(pz);
	    float EtaItem= E2eta(pz);
	    float PhiItem= E2phi(pz);

	    float RadOld=R2_bis(_Eta,_Phi,EtaItem,PhiItem);
	    float RadNew=R2_bis(NewJet.Eta(),NewJet.Phi(),EtaItem,PhiItem);
	    if (RadNew > RadOld) { 
	      it = NewJet._LItems.erase(it);
	    }
	    else {
	      _LItems.erase(where);
	      ++it;
	    }
	  }
	  else ++it;
	}
	return true;
      }
    }
  

    float dist_R2(TemporaryJet &jet) const {
      float deta= _Eta-jet.Eta();
      //float dphi= kinem::delta_phi(_Phi,jet.Phi());
      float dphi= inline_maths::delta_phi(_Phi,jet.Phi());
      return (deta*deta+dphi*dphi); 
    }
 
    bool ItemInJet(const CalItem* tw) const {
      typename list<const CalItem*>::const_iterator
	where= find(_LItems.begin(),_LItems.end(),tw);
      if(where != _LItems.end()) return true;
      else return false;
    }
  
    bool updateEtaPhiEt() { 
      float ETsum = 0.0;
      float ETAsum= 0.0;
      float PHIsum= 0.0;
      float Esum= 0.0;
      typename list<const CalItem*>::iterator it;
      for(it=_LItems.begin(); it!=_LItems.end(); it++) {
        float ETk = (*it)->pT();
        // now done in CalCell/CalTower if((*it)->E() < 0.0) ETk= -ETk;

        //float ETAk= (*it)->eta();
        //float PHIk= (*it)->phi();
        float pz[4];
        (*it)->p4vec(pz);
        float ETAk= E2eta(pz);
	// take care of the phi=0=2pi problem 
        float PHIk= E2phi(pz);
	//if(fabs(PHIk-_Phi) > kinem::TWOPI-fabs(PHIk-_Phi))
	if(fabs(PHIk-_Phi) > inline_maths::TWOPI-fabs(PHIk-_Phi))
	  {
	  if(_Phi < PHIk) 
	    {
	      //PHIk -= kinem::TWOPI;
	      PHIk -= inline_maths::TWOPI;
	    }
	  else
	    {
	      //PHIk += kinem::TWOPI;
	      PHIk += inline_maths::TWOPI;
	    }
	  }
	ETAsum+= ETAk*ETk;
	PHIsum+= PHIk*ETk;
	ETsum += ETk;
        // Esum+=(*it)->E(); use 4-momentum components 
	Esum+= pz[3];
      }
      if(ETsum <= 0.0) {
	_Eta= 0.0;
	_Phi= 0.0;
	_Et = 0.0;
        _E=0.;
	return false;
      }
      else {
         _Eta= ETAsum/ETsum;
         _Phi= PHIsum/ETsum; 
	 //if ( _Phi<0 ) _Phi+=kinem::TWOPI;
	 if ( _Phi<0 ) _Phi+=inline_maths::TWOPI;
         _Et = ETsum;
         _E  = Esum;
         return true;  
      }
    }

    void D0_Angle_updateEtaPhi() {
      float EXsum = 0.0;
      float EYsum = 0.0;
      float EZsum = 0.0;
      typename list<const CalItem*>::iterator it;
      for(it=_LItems.begin(); it!=_LItems.end(); it++) {
	float p[4];
	(*it)->p4vec(p);
	EXsum += p[0];
	EYsum += p[1];
	EZsum += p[2];
      }
      //_Phi=kinem::phi(EYsum,EXsum);
      _Phi=inline_maths::phi(EYsum,EXsum);
      //_Eta=kinem::eta(EXsum,EYsum,EZsum);
      _Eta=inline_maths::eta(EXsum,EYsum,EZsum);
    } 

    void getItems(list<const CalItem*> &ecv) const {
      ecv.clear(); //ls 27/Feb/2010
      typename list<const CalItem*>::const_iterator it;
      for(it=_LItems.begin(); it!=_LItems.end(); it++) {
	ecv.push_back(*it);
      }
    }

    float Eta() {return _Eta;}
    float Phi() {return _Phi;}
    float Et()  {return _Et;}
    float E()  {return _E;}
    list<const CalItem*> &LItems() {return _LItems;}
    
  private:
    list<const CalItem*> _LItems;
    float _Eta;
    float _Phi;
    float _Et;
    float _E;
  }; //class TemporaryJet

  void getItemsInCone(list<const CalItem*> &tlist, float etaJet, float phiJet,
  		      float cone_radius, float zvertex_in) const; 
  void getItemsInCone_bis(list<const CalItem*> &tlist, float etaJet, 
               float phiJet,float cone_radius, float zvertex_in) const; 
  
public:
  
  vector< TemporaryJet > TempColl;  

};
  /////////////////////////////////////////////////////////

//template < class CalItem,class CalItemAddress,class CalIClusterChunk >
template < class CalItem >
//void ConeClusterAlgo <CalItem,CalItemAddress,CalIClusterChunk >:: 
void ConeClusterAlgo <CalItem >:: 
getItemsInCone(list<const CalItem*> &tlist, float etaJet, float phiJet, 
	       float cone_radius, float zvertex_in) const {
//
// provide the list of Items (towers, Cells...) containing the energy from a 
// jet of a given cone size
//
  float ZVERTEX_MAX=200.;
  float DMIN=80.;
  float DMAX=360.;
  float THETA_margin=0.022;
   
  float zvertex=zvertex_in;
  float d1,d2;
  float phi_d1, phi_d2;
  float theta_E1, r1, r2, z1, z2;
  float theta_d1, theta_d2, eta_d1, eta_d2;

  // Ignore very large vertex positions
  if (fabs(zvertex) > ZVERTEX_MAX ) zvertex=0.0;
 
  if (zvertex >=0. ) {
    d1=fabs(DMIN-zvertex);
    d2=fabs(DMAX+zvertex);
  } else {
    d1=fabs(DMAX-zvertex);
    d2=fabs(DMIN+zvertex);
  }
  
  // calculate theta of physics cone and find which eta's this intercepts
  // a the maximum points
  phi_d1 = phiJet+cone_radius;
  //theta_E1 = kinem::theta(etaJet+cone_radius);
  theta_E1 = inline_maths::theta(etaJet+cone_radius);
  z1 = zvertex+d1*cos(theta_E1);
  r1 = d1*sin(theta_E1);

  phi_d2 = phiJet-cone_radius;
  //theta_E1 = kinem::theta(etaJet-cone_radius);
  theta_E1 = inline_maths::theta(etaJet-cone_radius);
  z2 = zvertex+d2*cos(theta_E1);
  r2 = d2*sin(theta_E1);

  // maximum spread in detector theta
  theta_d1 = atan2(r1, z1);
  theta_d2 = atan2(r2, z2);

  // make sure they stay in the calorimeter 
  theta_d1=max(theta_d1, THETA_margin);
  theta_d2=max(theta_d2, THETA_margin);
  //theta_d1=min(kinem::PI-(double)THETA_margin, (double)theta_d1);
  theta_d1=min(inline_maths::PI-(double)THETA_margin, (double)theta_d1);
  //theta_d2=min(kinem::PI-(double)THETA_margin, (double)theta_d2);
  theta_d2=min(inline_maths::PI-(double)THETA_margin, (double)theta_d2);

  //eta_d1 = kinem::eta(theta_d1);
  eta_d1 = inline_maths::eta(theta_d1);
  //eta_d2 = kinem::eta(theta_d2);
  eta_d2 = inline_maths::eta(theta_d2);


  typename list<const CalItem*>::iterator it;
  for (it=tlist.begin() ; it != tlist.end() ; ) {
    //float eta_cur= (*it)->eta();
    //float phi_cur= (*it)->phi();
    float pz[4];
    (*it)->p4vec(pz);
    float eta_cur= E2eta(pz);
    float phi_cur= E2phi(pz);

    bool accepted = eta_cur < eta_d1 && eta_cur > eta_d2;
    //if ( phi_d2>0 && phi_d1<kinem::TWOPI ) {
    if ( phi_d2>0 && phi_d1<inline_maths::TWOPI ) {
      accepted = accepted && phi_cur<phi_d1 && phi_cur>phi_d2;
    }
    else{ // case the cone overlap the phi=0=2pi line
      if ( phi_d2>0 ){
	accepted = accepted && 
	  //((phi_cur>phi_d2 && phi_cur<kinem::TWOPI) || phi_cur<phi_d1-kinem::TWOPI);
	  ((phi_cur>phi_d2 && phi_cur<inline_maths::TWOPI) || phi_cur<phi_d1-inline_maths::TWOPI);
      }
      else{
	accepted = accepted && 
	  //((phi_cur<phi_d1 && phi_cur>0) || phi_cur>phi_d2+kinem::TWOPI);
	  ((phi_cur<phi_d1 && phi_cur>0) || phi_cur>phi_d2+inline_maths::TWOPI);
      }
    }
    if ( ! accepted ) it = tlist.erase(it);
    else ++it;

  }
}
  /////////////////////////////////////////////////////////
//template < class CalItem,class CalItemAddress,class CalIClusterChunk >
template < class CalItem >
//void ConeClusterAlgo <CalItem,CalItemAddress,CalIClusterChunk >:: 
void ConeClusterAlgo <CalItem>:: 
getItemsInCone_bis(list<const CalItem*> &tlist, float etaJet, float phiJet, 
	       float cone_radius, float zvertex_in) const {
//
// provide the list of Items (towers, Cells...) containing the energy from a 
// jet of a given cone size
//
// WARNING: this is only to be used to compare to RUN I cone jets
// 
  float ZVERTEX_MAX=200.;
  float DMIN=80.;
  float DMAX=360.;
  float THETA_margin=0.022;
   
  float zvertex=zvertex_in;
  float d1,d2;
  float phi_d1, phi_d2;
  float theta_E1, r1, r2, z1, z2;
  float theta_d1, theta_d2, eta_d1, eta_d2;

  // Ignore very large vertex positions
  if (fabs(zvertex) > ZVERTEX_MAX ) zvertex=0.0;
 
  if (zvertex >=0. ) {
    d1=fabs(DMIN-zvertex);
    d2=fabs(DMAX+zvertex);
  } else {
    d1=fabs(DMAX-zvertex);
    d2=fabs(DMIN+zvertex);
  }
  
  // calculate theta of physics cone and find which eta's this intercepts
  // a the maximum points
  
  phi_d1 = phiJet+cone_radius;
  //theta_E1 = kinem::theta(etaJet+cone_radius);
  theta_E1 = inline_maths::theta(etaJet+cone_radius);
  z1 = zvertex+d1*cos(theta_E1);
  r1 = d1*sin(theta_E1);

  phi_d2 = phiJet-cone_radius;
  //theta_E1 = kinem::theta(etaJet-cone_radius);
  theta_E1 = inline_maths::theta(etaJet-cone_radius);
  z2 = zvertex+d2*cos(theta_E1);
  r2 = d2*sin(theta_E1);

  // maximum spread in detector theta

  theta_d1 = atan2(r1, z1);
  theta_d2 = atan2(r2, z2);

  // make sure they stay in the calorimeter 

  theta_d1=max(theta_d1, THETA_margin);
  theta_d2=max(theta_d2, THETA_margin);
  //theta_d1=min(kinem::PI-(double)THETA_margin, (double)theta_d1);
  theta_d1=min(inline_maths::PI-(double)THETA_margin, (double)theta_d1);
  //theta_d2=min(kinem::PI-(double)THETA_margin, (double)theta_d2);
  theta_d2=min(inline_maths::PI-(double)THETA_margin, (double)theta_d2);


  //eta_d1 = kinem::eta(theta_d1);
  eta_d1 = inline_maths::eta(theta_d1);
  //eta_d2 = kinem::eta(theta_d2);
  eta_d2 = inline_maths::eta(theta_d2);

  float signe;
 
  if( eta_d1>=0.0 ) signe= 1.0;
  else signe= -1.0;
  int ietaMAX= eta_d1/0.1+signe;
  if(fabs(eta_d1)>=4.45) ietaMAX= 37*signe; 
  else if(fabs(eta_d1)>=4.1) ietaMAX= 36*signe; 
  else if(fabs(eta_d1)>=3.7) ietaMAX= 35*signe; 
  else if(fabs(eta_d1)>=3.42) ietaMAX= 34*signe; 
  else if(fabs(eta_d1)>=3.2) ietaMAX= 33*signe; 
  
  if( eta_d2>=0.0 ) signe= 1.0;
  else signe= -1.0;
  int ietaMIN= eta_d2/0.1+signe;
  if(fabs(eta_d2)>=4.45) ietaMIN= 37*signe; 
  else if(fabs(eta_d2)>=4.1) ietaMIN= 36*signe; 
  else if(fabs(eta_d2)>=3.7) ietaMIN= 35*signe; 
  else if(fabs(eta_d2)>=3.42) ietaMIN= 34*signe; 
  else if(fabs(eta_d2)>=3.2) ietaMIN= 33*signe; 

  //int iphiMAX= 64*phi_d1/(2.*kinem::PI)+1;
  int iphiMAX= 64*phi_d1/(2.*inline_maths::PI)+1;
  //int iphiMIN= 64*phi_d2/(2.*kinem::PI)+1;
  int iphiMIN= 64*phi_d2/(2.*inline_maths::PI)+1;
 
  typename list<const CalItem*>::iterator it;
  for (it=tlist.begin() ; it != tlist.end() ; ) {
    //float eta_cur= (*it)->eta();
    //float phi_cur= (*it)->phi();
    int ieta= (*it)->address().ieta();
    int iphi= (*it)->address().iphi();
    
    bool accepted = ieta<ietaMAX && ieta>ietaMIN;
    if ( iphiMIN>0 && iphiMAX<=64 ) {
      accepted = accepted && iphi<iphiMAX && iphi>iphiMIN;
    }
    else{ // case the cone overlap the phi=0=2pi line
      if ( iphiMIN>0 ){
	accepted = accepted && 
	  ((iphi>iphiMIN && iphi<=64) || iphi<iphiMAX-64);
      }
      else{
	accepted = accepted && 
	  ((iphi<iphiMAX && iphi>0) || iphi>iphiMIN+64);
      }
    }
    if ( ! accepted ) it = tlist.erase(it);
    else ++it;
    
  }
}
  /////////////////////////////////////////////////////////
//template < class CalItem,class CalItemAddress,class CalIClusterChunk >
template < class CalItem >
//inline void ConeClusterAlgo <CalItem,CalItemAddress,CalIClusterChunk >:: 
inline void ConeClusterAlgo <CalItem >:: 
print(ostream &os) const {
    os<<endl<<" CONE ALGORITHM, cone radius= "<<_CONErad<<endl<<
    " min E_T fraction= "<<_JETmne<<endl<<
    " minimum Delta_R separation between cones= "<<_TWOrad<<endl<<
    " shared E_T fraction threshold for combining jets= "<<_SPLifr<<endl;
}
  /////////////////////////////////////////////////////////

//template < class CalItem,class CalItemAddress,class CalIClusterChunk >
template < class CalItem >
//void ConeClusterAlgo <CalItem,CalItemAddress,CalIClusterChunk >:: 
void ConeClusterAlgo <CalItem >:: 
makeClusters(//const EnergyClusterReco* r,
	     std::list<CalItem> &jets,
	     list<const CalItem*> &itemlist, float Zvertex 
	     //, const EnergyClusterCollection<CalItemAddress> &preclu,
	     //CalIClusterChunk* chunkptr
	     //) const {
	     ) {

  // create an energy cluster collection for jets 
  //EnergyClusterCollection<CalItemAddress>* ptrcol;
  //r->createClusterCollection(chunkptr, ptrcol);
  std::vector<const CalItem*> ecv;
  for ( typename std::list<const CalItem*>::iterator it = itemlist.begin(); 
        it != itemlist.end(); it++) {
    ecv.push_back(*it);
  }


  // Initialize
  float Rcut= 1.E-06;
  if(_Increase_Delta_R) Rcut= 1.E-04;
  bool nojets;

  //vector< TemporaryJet > TempColl;  
  list< pair<float,float> > LTrack;

  // get a vector with pointers to EnergyCluster in the collection
  //vector<const EnergyCluster<CalItemAddress>*> ecv;
  //preclu.getClusters(ecv);

  // loop over all preclusters
  //typename vector<const EnergyCluster<CalItemAddress>*>::iterator jclu;
  typename std::vector<const CalItem*>::iterator jclu;
  for( jclu=ecv.begin(); jclu!=ecv.end(); jclu++ ) {
    ////const EnergyCluster<CalItemAddress>* ptr= *jclu;
    const CalItem* ptr= *jclu;
    //float PHIst= ptr->phi();
    //float ETAst= ptr->eta();
    float pz[4];
    ptr->p4vec(pz);
    float ETAst= E2eta(pz);
    float PHIst= E2phi(pz);

    //cout << "seed 4-vec ";
    //for ( int i = 0; i < 4; i++) cout << pz[i] << " ";
    //cout << endl;

    nojets= false;
    // check to see if precluster is too close to a found jet
    if(_Kill_Far_Clusters) {
      list< pair<float,float> >::iterator kj;
      for(kj=LTrack.begin(); kj!=LTrack.end(); kj++) {
	if(DELTA_r((*kj).first,ETAst,(*kj).second,PHIst)<_Far_Ratio*_CONErad) {
	  nojets= true;
	  //cout << "seed too close ! skip " << endl;
	  break;
	}
      }
    }
    if( nojets==false ) {
      TemporaryJet TJet(ETAst,PHIst);
      list< pair<int,float> > JETshare;

      // start of cone building loop
      int trial= 0;
      do{  
        trial++;
	//cout << " trial " << trial << endl;

	ETAst= TJet.Eta();
        PHIst= TJet.Phi();
        TJet.erase();

	//if(PHIst > kinem::TWOPI) PHIst= PHIst-kinem::TWOPI;
	if(PHIst > inline_maths::TWOPI) PHIst= PHIst-inline_maths::TWOPI;
	//if(PHIst < 0.0  ) PHIst= PHIst+kinem::TWOPI;
	if(PHIst < 0.0  ) PHIst= PHIst+inline_maths::TWOPI;
	//if( PHIst>kinem::TWOPI || PHIst<0.0 ) {
	if( PHIst>inline_maths::TWOPI || PHIst<0.0 ) {
          TJet.setEtaPhiEt(0.0,0.0,0.0);
	  break; // end loop do (illegal jet PHI) 
        }
	TJet.setEtaPhiEt(ETAst,PHIst,0.0);

	// calculate eta & phi limits for cone
        list<const CalItem*> Twlist(itemlist);

	getItemsInCone(Twlist,ETAst,PHIst,_CONErad,Zvertex); 
	//  only to compare with RUN I cone jets !   getItemsInCone_bis(Twlist,ETAst,PHIst,_CONErad,Zvertex); 

	// loop over the possible items for this cone
        typename list<const CalItem*>::iterator tk;
        for( tk= Twlist.begin(); tk!=Twlist.end(); tk++ ) {
	  float ETk =(*tk)->pT();
          // now done in CalCell/CalTower if((*tk)->E() < 0.0) ETk= -ETk;

          if( ETk > _Eitem_Negdrop ) { 
	    //float ETAk=(*tk)->eta();
	    //float PHIk=(*tk)->phi();
            float pz[4];
            (*tk)->p4vec(pz);
            float ETAk= E2eta(pz);
            float PHIk= E2phi(pz);

	    float dphi= fabs(PHIk-PHIst);
	    //if(dphi > kinem::TWOPI-dphi) {
	    if(dphi > inline_maths::TWOPI-dphi) {
              //if(PHIst < PHIk) PHIk= PHIk-kinem::TWOPI;
              if(PHIst < PHIk) PHIk= PHIk-inline_maths::TWOPI;
              //else PHIk= PHIk+kinem::TWOPI; 
              else PHIk= PHIk+inline_maths::TWOPI; 
            }

	    if( R2_bis(ETAk,PHIk,ETAst,PHIst) <= _CONErad*_CONErad ) { 
	      TJet.addItem(*tk);
	    }
	  }
	}// end loop tk
 
	if(TJet.updateEtaPhiEt()==false) {
	  //cout << " negative E jet ! drop " << endl;
	  break;
	}

	// require some minimum ET on every iteration
	if(_Jet_Et_Min_On_Iter) {
	  if( TJet.Et() < _JETmne*_Et_Min_Ratio ) {
	    //cout << " too low ET jet ! drop " << endl;
	    break; // end loop trial
          }
        }
  
	//cout << " R2 = " << R2_bis(TJet.Eta(),TJet.Phi(),ETAst,PHIst) << 
	//  " Rcut = " << Rcut << endl;
      }while(R2_bis(TJet.Eta(),TJet.Phi(),ETAst,PHIst)>=Rcut && trial<=50);
 
      if( TJet.Et() >= _JETmne ) {
	//cout << " jet accepted will check for overlaps " << endl; 
	if(_D0_Angle) TJet.D0_Angle_updateEtaPhi();
	//cout << "  after TJet.D0_Angle_updateEtaPhi() " << endl;
	
        // item also in another jet
        list<const CalItem*> Lst;
        TJet.getItems(Lst);
        typename list<const CalItem*>::iterator tk;
        for(tk=Lst.begin(); tk!=Lst.end(); tk++) {
	  float ETk=(*tk)->pT();
          // now done in CalCell/CalTower if((*tk)->E() < 0.0) ETk= -ETk;
          for(unsigned int kj=0; kj<TempColl.size(); kj++) {
	    if(TempColl[kj].ItemInJet(*tk)==true) {
	      list< pair<int,float> >::iterator pit;
              bool jetok= false;
              for(pit=JETshare.begin(); pit!=JETshare.end();pit++) {
                if((*pit).first == (int) kj) {
                  jetok= true;
                  (*pit).second+= ETk;
                  break;
                }
              }
              if(jetok==false) JETshare.push_back(make_pair(kj,ETk));
            }
          }
        }
	
	if(JETshare.size() >0) {
	  list< pair<int,float> >::iterator pit;
	  float Ssum= 0.0;
	  list< pair<int,float> >::iterator pitMAX=JETshare.begin();
	  for(pit=JETshare.begin(); pit!=JETshare.end(); pit++) {
	    Ssum+= (*pit).second;
	    if((*pit).second > (*pitMAX).second) pitMAX= pit;
	  }

          //int IJET= (*pitMAX).first;
	  bool splshr;
	  float Eleft= fabs(TJet.Et()-Ssum);
	  float Djets= TempColl[(*pitMAX).first].dist_R2(TJet);
	  if(Djets <= _TWOrad || Eleft <= _Thresh_Diff_Et) { 
	    TJet.erase();
	    splshr= false;
	  }
	  else {
	    float SharedFr=Ssum/min(TempColl[(*pitMAX).first].Et(),TJet.Et());
	    if(JETshare.size() >1) {
	      typename list<const CalItem*>::iterator tk;
	      for(tk=TJet.LItems().begin(); tk!=TJet.LItems().end(); ) {
                bool found = false;
		list< pair<int,float> >::iterator pit;
		for(pit=JETshare.begin(); pit!=JETshare.end();pit++) {
		  if((*pit).first!=(*pitMAX).first) { 
		    if(TempColl[(*pit).first].ItemInJet(*tk)==true) {
		      tk = TJet.LItems().erase(tk);
		      found = true;
		      break;
		    }
		  }
		}
		if ( !found ) ++tk;
	      }
	    }

	    splshr= TempColl[(*pitMAX).first].share_jets(TJet,SharedFr,_SPLifr);

	    if( splshr==true ) {
	      //cout << " jet splitted due to overlaps " << endl;
	      TempColl[(*pitMAX).first].updateEtaPhiEt();
	      TJet.updateEtaPhiEt();
	      if(_D0_Angle) TJet.D0_Angle_updateEtaPhi();
	      if(_D0_Angle) TempColl[(*pitMAX).first].D0_Angle_updateEtaPhi();
	      TempColl.push_back(TJet);  
	      LTrack.push_back(make_pair(TJet.Eta(),TJet.Phi()));
	    }
	    else {
	      //cout << " jet merged due to overlaps " << endl;
	      TempColl[(*pitMAX).first].updateEtaPhiEt();
	      if(_D0_Angle) TempColl[(*pitMAX).first].D0_Angle_updateEtaPhi();
	    }  
	  }
	}
	else {
	  TJet.updateEtaPhiEt();
	  if(_D0_Angle) TJet.D0_Angle_updateEtaPhi();
          TempColl.push_back(TJet);  
          LTrack.push_back(make_pair(TJet.Eta(),TJet.Phi()));
	}
      } //JETmne
    } //nojets
  }// end loop jclu 

  for(unsigned int i=0; i<TempColl.size(); i++) {
    //EnergyCluster<CalItemAddress>* ptrclu;
    CalItem ptrclu;
    //r->createCluster(ptrcol,ptrclu);
    list<const CalItem*> Vi;
    TempColl[i].getItems(Vi);
    typename list<const CalItem*>::iterator it;
    for(it=Vi.begin(); it!=Vi.end(); it++) {
      const CalItem* ptr= *it;
      //CalItemAddress addr= ptr->address();
      float p[4];
      ptr->p4vec(p);
      //float emE= ptr->emE();
      //r->addClusterItem(ptrclu,addr,p,emE);
      ptrclu.Add(*ptr);
    }
    jets.push_back(ptrclu);
  }

}// end 

} //namespace d0runi

FASTJET_END_NAMESPACE

#endif  //  CONECLUSTERALGO_H



