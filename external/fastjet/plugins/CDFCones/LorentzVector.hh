#ifndef _LORENTZ_VECTOR_HH_
#define _LORENTZ_VECTOR_HH_

//----------------------------------------------------------------------
// This file distributed with FastJet has been obtained from
// http://www.pa.msu.edu/~huston/Les_Houches_2005/JetClu+Midpoint-StandAlone.tgz
//
// Permission to distribute it with FastJet has been granted by Joey
// Huston (see the COPYING file in the main FastJet directory for
// details).
// Changes from the original file are listed below.
//----------------------------------------------------------------------

// History of changes compared to the original LorentzVector.hh file
// 
// 2011-11-14  Gregory Soyez  <soyez@fastjet.fr>
//
//        * removed some harmless warnings coming with the -Wshadow gcc option
//
// 2009-01-17  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * put the code in the fastjet::cdf namespace
// 
// 2008-01-15  Gregory Soyez  <soyez@fastjet.fr>
// 
// 	  * fixed issues with compilation under VC (definition of M_PI)
// 
// 2007-02-21  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * added option of choosing the scale used in the split-merge
//          procedure (pt [default], Et or mt)
// 
// 2006-09-24  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * added JetClu+MidPoint to FastJet

#include <cmath>

#ifndef M_PI
#define M_PI  3.141592653589793238462643383279502884197 
#endif

#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace cdf{

class LorentzVector
{
 public:

  double px,py,pz,E;

  LorentzVector(): px(0), py(0), pz(0), E(0) {}
  LorentzVector(double p1, double p2, double p3, double p0): px(p1), py(p2), pz(p3), E(p0) {}
  LorentzVector(const LorentzVector& lv): px(lv.px), py(lv.py), pz(lv.pz), E(lv.E) {}
  double p()   const {return sqrt(px*px + py*py + pz*pz);}
  double pt()  const {return sqrt(px*px + py*py);}
  double mt()  const {return sqrt((E-pz)*(E+pz));}
  double y()   const {return 0.5*log((E + pz)/(E - pz));}
  double Et()  const {return E/p()*pt();}
  double eta() const {return 0.5*log((p() + pz)/(p() - pz));}
  double phi() const
  {
    double r = atan2(py,px);
    if(r < 0)
      r += 2*M_PI;
    return r;
  }
  void add(LorentzVector v)
  {
    px += v.px;
    py += v.py;
    pz += v.pz;
    E  += v.E;
  }
  void subtract(LorentzVector v)
  {
    px -= v.px;
    py -= v.py;
    pz -= v.pz;
    E  -= v.E;
  }
  bool isEqual(LorentzVector v)
  {
    return px == v.px && py == v.py && pz == v.pz && E == v.E;
  }
};

}  // namespace cdf

FASTJET_END_NAMESPACE

#endif
