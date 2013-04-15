#ifndef  D0RunIIconeJets_HepEntity_class
#define  D0RunIIconeJets_HepEntity_class

#include "inline_maths.h"

#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace d0{

//Author: Lars Sonnenschein 28/Mar/2007
//This is an example class fulfilling the minimal requirements needed by the
//D0 RunII cone jet algorithm implementation, which is an inlined template class
//
//
// This file is distributed with FastJet under the terms of the GNU
// General Public License (v2). Permission to do so has been granted
// by Lars Sonnenschein and the D0 collaboration (see COPYING for
// details)
//
// History of changes in FastJet compared tothe original version of
// HepEntity.h
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
//        * added an index member

class HepEntity {

 public:

  HepEntity() {
    E=0.;
    px=0.;
    py=0.;
    pz=0.;
    index = -1;
    return;
  }


  HepEntity(double E_in, double px_in, double py_in, double pz_in, int index_in = -1) : 
    E(E_in), px(px_in), py(py_in), pz(pz_in), index(index_in) {
    return;
  }


    HepEntity(const HepEntity& in) : E(in.E), px(in.px), py(in.py), pz(in.pz), index(in.index) {
    return;
  }

  
  inline double y() const {
    return inline_maths::y(E,pz);
  }


  inline double phi() const {
     return inline_maths::phi(px,py);
  }


  inline double pT() const {
     return sqrt(inline_maths::sqr(px)+inline_maths::sqr(py));
  }


  inline void p4vec(float* p) const {
    p[0] = px;
    p[1] = py;
    p[2] = pz;
    p[3] = E;
    return;
  }

  inline void Add(const HepEntity el) {
    E += el.E;
    px += el.px;
    py += el.py;
    pz += el.pz;
    return;
  }

  inline void Fill(double E_in, double px_in, double py_in, double pz_in, int index_in = -1) {
    E = E_in;
    px = px_in;
    py = py_in;
    pz = pz_in;
    index = index_in;
    return;
  }


  double E;
  double px;
  double py;
  double pz;
  int    index;

 private:



};
//end of class HepEntity;


}  // namespace d0

FASTJET_END_NAMESPACE

#endif
