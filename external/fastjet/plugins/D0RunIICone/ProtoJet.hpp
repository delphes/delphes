#ifndef  D0RunIIconeJets_PROTOJET
#define  D0RunIIconeJets_PROTOJET
// ---------------------------------------------------------------------------
// ProtoJet.hpp
//
// Created: 28-JUL-2000 Francois Touze (+ Laurent Duflot)
//
// Purpose: Implements a proto-jet object that is used as input by the 
//   Improved Legacy Cone Algorithm split/merge algo.
//
// Modified:
//    9-Aug-2000  Laurent Duflot
//     + save the initial stable cone ET before split/merge
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
// ProtoJet.hpp
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
//        * replaced make_pair by std::make_pair
//
// ---------------------------------------------------------------------------
 
//#include "kinem_util/AnglesUtil.hpp"
//#include "energycluster/ConeJetInfo.hpp"
#include "ConeJetInfo.hpp"
#include <list>
#include <cmath>

#include "inline_maths.h" //ls

#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace d0{

using namespace inline_maths;
using namespace D0RunIIconeJets_CONEJETINFO;


inline float RD2(float y1,float phi1,float y2,float phi2) 
{
  float dphi= delta_phi(phi1,phi2);
  return (y1-y2)*(y1-y2)+dphi*dphi; 
}

inline float RDelta(float y1,float phi1,float y2,float phi2) 
{
  float dphi= delta_phi(phi1,phi2);
  return sqrt((y1-y2)*(y1-y2)+dphi*dphi); 
}

inline float P2y(float* p4vec) {
  return y(p4vec[3],p4vec[2]);
}

inline float P2phi(float* p4vec) {
  return phi(p4vec[0],p4vec[1]);
}

///////////////////////////////////////////////////////////////////////////////
template <class Item>
class ProtoJet {

public :

  ProtoJet(float seedET);
  ProtoJet(float seedET,float y,float phi);
  ProtoJet(const ProtoJet<Item>& pj);
  ~ProtoJet() {;}

  void addItem(const Item* tw); 
  void setJet(float y,float phi,float pT); 
  void updateJet();
  void erase();

  float y() const; 
  float phi() const;
  float pT() const;
  const ConeJetInfo & info() const;
  const std::list<const Item*>& LItems() const;

  void print(std::ostream &os) const;

  // actions to be taken when the jet is a stable cone
  void NowStable();
  // declare the jet to have been splitted
  void splitted(){_info.splitted();};
  // declare the jet to have been merged
  void merged(){_info.merged();};
protected :

  std::list<const Item*> _LItems;
  float _y;
  float _phi;
  float _pT;
  ConeJetInfo _info;

};
///////////////////////////////////////////////////////////////////////////////
template<class Item>
ProtoJet<Item>::ProtoJet(float seedET) : _LItems(), _info(seedET) {
    _y  = 0.0;
    _phi= 0.0;
    _pT = 0.0;
}

template<class Item>
ProtoJet<Item>::ProtoJet(float seedET,float y_in,float phi_in) :  _LItems(), _info(seedET) { 
  _y  = y_in; 
  _phi= phi_in;
  _pT = 0.0;
}

template<class Item>
ProtoJet<Item>::ProtoJet(const ProtoJet<Item>& pj): _y(pj._y), 
						    _phi(pj._phi), _pT(pj._pT),
                                                    _info(pj._info)
{ 
  typename std::list<const Item*>::const_iterator it;
  for(it = pj._LItems.begin(); it != pj._LItems.end(); ++it) { 
    _LItems.push_back(*it);
  }
}

template<class Item>
void ProtoJet<Item>::addItem(const Item* tw) {
  _LItems.push_back(tw);
}

template<class Item>
void ProtoJet<Item>::setJet(float y_in,float phi_in,float pT_in) {
  _y  = y_in;
  _phi= phi_in;
  _pT = pT_in;
}

template<class Item>
void ProtoJet<Item>::updateJet() { 
  //float ETsum = 0.0;
  //float ysum  = 0.0;
  //float PHIsum= 0.0;
  float p[4] = {0.,0.,0.,0.};
  typename std::list<const Item*>::iterator it;
  for(it = _LItems.begin(); it != _LItems.end(); ++it) 
  {
    float pk[4];
    (*it)->p4vec(pk);
    //cout << "updateJet: px=" << pk[0] << " py=" << pk[1] << " pz=" << pk[2] << " E=" << pk[3] << endl; 
    for ( int i = 0; i < 4 ; ++i) p[i] += pk[i];
  }
  _y = P2y(p);
  _phi = P2phi(p);
  _pT = sqrt(p[0]*p[0] + p[1]*p[1]);
  if ( p[3] < 0. ) _pT = - _pT;

}

template<class Item>
void ProtoJet<Item>::erase() {
  _LItems.erase(_LItems.begin(),_LItems.end());
  _y  = 0.0;
  _phi= 0.0;
  _pT = 0.0; 
  // _info is not modified in order to keep split/merge history
}

// actions to be taken when the jet is a stable cone
template<class Item>
void ProtoJet<Item>::NowStable() {
  _info.initialET(_pT);
}

template<class Item>
void ProtoJet<Item>::print(std::ostream& os) const {
  os<<"y phi Et = ("<<_y<<", "<<_phi<<", "<<this->_Et<<")"<<std::endl;
  os<< " members= " << std::endl;
  typename std::list<const Item*>::const_iterator i;
  for(i = _LItems.begin(); i != _LItems.end(); ++i)
    (*i)->print(os);
  os << std::endl;
}

template<class Item>
inline float ProtoJet<Item>::y() const{
  return _y;
}

template<class Item>
inline float ProtoJet<Item>::phi() const{
  return _phi;
}

template<class Item>
inline float ProtoJet<Item>::pT() const{
  return _pT;
}
template<class Item>
inline const ConeJetInfo & ProtoJet<Item>::info() const{
  return _info;
}

template<class Item>   
inline const std::list<const Item*>& ProtoJet<Item>::LItems() const{
  return _LItems;
}
///////////////////////////////////////////////////////////////////////////////

}  // namespace d0

FASTJET_END_NAMESPACE

#endif
