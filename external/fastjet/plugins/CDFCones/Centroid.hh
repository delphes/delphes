#ifndef _CENTROID_HH_
#define _CENTROID_HH_

//----------------------------------------------------------------------
// This file distributed with FastJet has been obtained from
// http://www.pa.msu.edu/~huston/Les_Houches_2005/JetClu+Midpoint-StandAlone.tgz
//
// Permission to distribute it with FastJet has been granted by Joey
// Huston (see the COPYING file in the main FastJet directory for
// details).
// Changes from the original file are listed below.
//----------------------------------------------------------------------

// History of changes compared to the original Centroid.hh file
// 
// 2009-01-17  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * put the code in the fastjet::cdf namespace
// 
// 2008-01-15  Gregory Soyez  <soyez@fastjet.fr>
// 
// 	  * fixed issues with compilation under VC (definition of M_PI)
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

class Centroid
{
 public:

  double Et,eta,phi;

  Centroid(): Et(0), eta(0), phi(0) {}
  Centroid(double centroidEt, double centroidEta, double centroidPhi): Et(centroidEt), eta(centroidEta), phi(centroidPhi) {}
  Centroid(const Centroid& c): Et(c.Et), eta(c.eta), phi(c.phi) {}
  void add(Centroid c)
  {
    double newEt = Et + c.Et;
    eta = (Et*eta + c.Et*c.eta)/newEt;
    double dPhi = c.phi - phi;
    if(dPhi > M_PI)
      dPhi -= 2*M_PI;
    else if(dPhi < -M_PI)
      dPhi += 2*M_PI;
    phi += dPhi*c.Et/newEt;
    while(phi < 0)
      phi += 2*M_PI;
    while(phi >= 2*M_PI)
      phi -= 2*M_PI;
    Et = newEt;
  }
  void subtract(Centroid c)
  {
    double newEt = Et - c.Et;
    eta = (Et*eta - c.Et*c.eta)/newEt;
    double dPhi = c.phi - phi;
    if(dPhi > M_PI)
      dPhi -= 2*M_PI;
    else if(dPhi < -M_PI)
      dPhi += 2*M_PI;
    phi -= dPhi*c.Et/newEt;
    while(phi < 0)
      phi += 2*M_PI;
    while(phi >= 2*M_PI)
      phi -= 2*M_PI;
    Et = newEt;
  }
  bool isEqual(Centroid c)
  {
    return Et == c.Et && eta == c.eta && phi == c.phi;
  }
};

}  // namespace cdf

FASTJET_END_NAMESPACE

#endif
