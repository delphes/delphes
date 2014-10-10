//FJSTARTHEADER
// $Id: PseudoJet.cc 3652 2014-09-03 13:31:13Z salam $
//
// Copyright (c) 2005-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
//  FastJet as part of work towards a scientific publication, please
//  quote the version you use and include a citation to the manual and
//  optionally also to hep-ph/0512210.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//FJENDHEADER


#include "fastjet/Error.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#ifndef __FJCORE__
#include "fastjet/ClusterSequenceAreaBase.hh"
#endif  // __FJCORE__
#include "fastjet/CompositeJetStructure.hh"
#include<valarray>
#include<iostream>
#include<sstream>
#include<cmath>
#include<algorithm>
#include <cstdarg>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;


//----------------------------------------------------------------------
// another constructor...
PseudoJet::PseudoJet(const double px_in, const double py_in, const double pz_in, const double E_in) {
  
  _E  = E_in ;
  _px = px_in;
  _py = py_in;
  _pz = pz_in;

  this->_finish_init();

  // some default values for the history and user indices
  _reset_indices();

}


//----------------------------------------------------------------------
/// do standard end of initialisation
void PseudoJet::_finish_init () {
  _kt2 = this->px()*this->px() + this->py()*this->py();
  _phi = pseudojet_invalid_phi;
  // strictly speaking, _rap does not need initialising, because
  // it's never used as long as _phi == pseudojet_invalid_phi
  // (and gets set when _phi is requested). However ATLAS
  // 2013-03-28 complained that they sometimes have NaN's in
  // _rap and this interferes with some of their internal validation. 
  // So we initialise it; penalty is about 0.3ns per PseudoJet out of about
  // 10ns total initialisation time (on a intel Core i7 2.7GHz)
  _rap = pseudojet_invalid_rap;
}

//----------------------------------------------------------------------
void PseudoJet::_set_rap_phi() const {

  if (_kt2 == 0.0) {
    _phi = 0.0; } 
  else {
    _phi = atan2(this->py(),this->px());
  }
  if (_phi < 0.0) {_phi += twopi;}
  if (_phi >= twopi) {_phi -= twopi;} // can happen if phi=-|eps<1e-15|?
  if (this->E() == abs(this->pz()) && _kt2 == 0) {
    // Point has infinite rapidity -- convert that into a very large
    // number, but in such a way that different 0-pt momenta will have
    // different rapidities (so as to lift the degeneracy between
    // them) [this can be relevant at parton-level]
    double MaxRapHere = MaxRap + abs(this->pz());
    if (this->pz() >= 0.0) {_rap = MaxRapHere;} else {_rap = -MaxRapHere;}
  } else {
    // get the rapidity in a way that's modestly insensitive to roundoff
    // error when things pz,E are large (actually the best we can do without
    // explicit knowledge of mass)
    double effective_m2 = max(0.0,m2()); // force non tachyonic mass
    double E_plus_pz    = _E + abs(_pz); // the safer of p+, p-
    // p+/p- = (p+ p-) / (p-)^2 = (kt^2+m^2)/(p-)^2
    _rap = 0.5*log((_kt2 + effective_m2)/(E_plus_pz*E_plus_pz));
    if (_pz > 0) {_rap = - _rap;}
  }

}


//----------------------------------------------------------------------
// return a valarray four-momentum
valarray<double> PseudoJet::four_mom() const {
  valarray<double> mom(4);
  mom[0] = _px;
  mom[1] = _py;
  mom[2] = _pz;
  mom[3] = _E ;
  return mom;
}

//----------------------------------------------------------------------
// Return the component corresponding to the specified index.
// taken from CLHEP
double PseudoJet::operator () (int i) const {
  switch(i) {
  case X:
    return px();
  case Y:
    return py();
  case Z:
    return pz();
  case T:
    return e();
  default:
    ostringstream err;
    err << "PseudoJet subscripting: bad index (" << i << ")";
    throw Error(err.str());
  }
  return 0.;
}  

//----------------------------------------------------------------------
// return the pseudorapidity
double PseudoJet::pseudorapidity() const {
  if (px() == 0.0 && py() ==0.0) return MaxRap;
  if (pz() == 0.0) return 0.0;

  double theta = atan(perp()/pz());
  if (theta < 0) theta += pi;
  return -log(tan(theta/2));
}

//----------------------------------------------------------------------
// return "sum" of two pseudojets
PseudoJet operator+ (const PseudoJet & jet1, const PseudoJet & jet2) {
  //return PseudoJet(jet1.four_mom()+jet2.four_mom());
  return PseudoJet(jet1.px()+jet2.px(),
		   jet1.py()+jet2.py(),
		   jet1.pz()+jet2.pz(),
		   jet1.E() +jet2.E()  );
} 

//----------------------------------------------------------------------
// return difference of two pseudojets
PseudoJet operator- (const PseudoJet & jet1, const PseudoJet & jet2) {
  //return PseudoJet(jet1.four_mom()-jet2.four_mom());
  return PseudoJet(jet1.px()-jet2.px(),
		   jet1.py()-jet2.py(),
		   jet1.pz()-jet2.pz(),
		   jet1.E() -jet2.E()  );
} 

//----------------------------------------------------------------------
// return the product, coeff * jet
PseudoJet operator* (double coeff, const PseudoJet & jet) {
  // see the comment in operator*= about ensuring valid rap phi
  // before a multiplication to handle case of multiplication by
  // zero, while maintaining rapidity and phi
  jet._ensure_valid_rap_phi(); 
  //return PseudoJet(coeff*jet.four_mom());
  // the following code is hopefully more efficient
  PseudoJet coeff_times_jet(jet);
  coeff_times_jet *= coeff;
  return coeff_times_jet;
} 

//----------------------------------------------------------------------
// return the product, coeff * jet
PseudoJet operator* (const PseudoJet & jet, double coeff) {
  return coeff*jet;
} 

//----------------------------------------------------------------------
// return the ratio, jet / coeff
PseudoJet operator/ (const PseudoJet & jet, double coeff) {
  return (1.0/coeff)*jet;
} 

//----------------------------------------------------------------------
/// multiply the jet's momentum by the coefficient
void PseudoJet::operator*=(double coeff) {
  // operator*= aims to maintain the rapidity and azimuth
  // for the PseudoJet; if they have already been evaluated
  // this is fine, but if they haven't and coeff is sufficiently
  // small as to cause a zero or underflow result, then a subsequent
  // invocation of rap or phi will lead to a non-sensical result. 
  // So, here, we preemptively ensure that rapidity and phi
  // are correctly cached
  _ensure_valid_rap_phi(); 
  _px *= coeff;
  _py *= coeff;
  _pz *= coeff;
  _E  *= coeff;
  _kt2*= coeff*coeff;
  // phi and rap are unchanged
}

//----------------------------------------------------------------------
/// divide the jet's momentum by the coefficient
void PseudoJet::operator/=(double coeff) {
  (*this) *= 1.0/coeff;
}


//----------------------------------------------------------------------
/// add the other jet's momentum to this jet
void PseudoJet::operator+=(const PseudoJet & other_jet) {
  _px += other_jet._px;
  _py += other_jet._py;
  _pz += other_jet._pz;
  _E  += other_jet._E ;
  _finish_init(); // we need to recalculate phi,rap,kt2
}


//----------------------------------------------------------------------
/// subtract the other jet's momentum from this jet
void PseudoJet::operator-=(const PseudoJet & other_jet) {
  _px -= other_jet._px;
  _py -= other_jet._py;
  _pz -= other_jet._pz;
  _E  -= other_jet._E ;
  _finish_init(); // we need to recalculate phi,rap,kt2
}

//----------------------------------------------------------------------
bool operator==(const PseudoJet & a, const PseudoJet & b) {
  if (a.px() != b.px()) return false;
  if (a.py() != b.py()) return false;
  if (a.pz() != b.pz()) return false;
  if (a.E () != b.E ()) return false;
  
  if (a.user_index()    != b.user_index()) return false;
  if (a.cluster_hist_index() != b.cluster_hist_index()) return false;
  if (a.user_info_ptr() != b.user_info_ptr()) return false;
  if (a.structure_ptr() != b.structure_ptr()) return false;

  return true;
}

//----------------------------------------------------------------------
// check if the jet has zero momentum
bool operator==(const PseudoJet & jet, const double val) {
  if (val != 0) 
    throw Error("comparing a PseudoJet with a non-zero constant (double) is not allowed.");
  return (jet.px() == 0 && jet.py() == 0 && 
	  jet.pz() == 0 && jet.E() == 0);
}



//----------------------------------------------------------------------
/// transform this jet (given in the rest frame of prest) into a jet
/// in the lab frame 
//
// NB: code adapted from that in herwig f77 (checked how it worked
// long ago)
PseudoJet & PseudoJet::boost(const PseudoJet & prest) {
  
  if (prest.px() == 0.0 && prest.py() == 0.0 && prest.pz() == 0.0) 
    return *this;

  double m_local = prest.m();
  assert(m_local != 0);

  double pf4  = (  px()*prest.px() + py()*prest.py()
                 + pz()*prest.pz() + E()*prest.E() )/m_local;
  double fn   = (pf4 + E()) / (prest.E() + m_local);
  _px +=  fn*prest.px();
  _py +=  fn*prest.py();
  _pz +=  fn*prest.pz();
  _E = pf4;

  _finish_init(); // we need to recalculate phi,rap,kt2
  return *this;
}


//----------------------------------------------------------------------
/// transform this jet (given in lab) into a jet in the rest
/// frame of prest  
//
// NB: code adapted from that in herwig f77 (checked how it worked
// long ago)
PseudoJet & PseudoJet::unboost(const PseudoJet & prest) {
  
  if (prest.px() == 0.0 && prest.py() == 0.0 && prest.pz() == 0.0) 
    return *this;

  double m_local = prest.m();
  assert(m_local != 0);

  double pf4  = ( -px()*prest.px() - py()*prest.py()
                 - pz()*prest.pz() + E()*prest.E() )/m_local;
  double fn   = (pf4 + E()) / (prest.E() + m_local);
  _px -=  fn*prest.px();
  _py -=  fn*prest.py();
  _pz -=  fn*prest.pz();
  _E = pf4;

  _finish_init(); // we need to recalculate phi,rap,kt2
  return *this;
}


//----------------------------------------------------------------------
/// returns true if the momenta of the two input jets are identical
bool have_same_momentum(const PseudoJet & jeta, const PseudoJet & jetb) {
  return jeta.px() == jetb.px()
    &&   jeta.py() == jetb.py()
    &&   jeta.pz() == jetb.pz()
    &&   jeta.E()  == jetb.E();
}

//----------------------------------------------------------------------
void PseudoJet::set_cached_rap_phi(double rap_in, double phi_in) {
  _rap = rap_in; _phi = phi_in;
  if (_phi >= twopi) _phi -= twopi;
  if (_phi < 0)      _phi += twopi;
}

//----------------------------------------------------------------------
void PseudoJet::reset_momentum_PtYPhiM(double pt_in, double y_in, double phi_in, double m_in) {
  assert(phi_in < 2*twopi && phi_in > -twopi);
  double ptm = (m_in == 0) ? pt_in : sqrt(pt_in*pt_in+m_in*m_in);
  double exprap = exp(y_in);
  double pminus = ptm/exprap;
  double pplus  = ptm*exprap;
  double px_local = pt_in*cos(phi_in);
  double py_local = pt_in*sin(phi_in);
  reset_momentum(px_local,py_local,0.5*(pplus-pminus),0.5*(pplus+pminus));
  set_cached_rap_phi(y_in,phi_in);
}

//----------------------------------------------------------------------
/// return a pseudojet with the given pt, y, phi and mass
PseudoJet PtYPhiM(double pt, double y, double phi, double m) {
  assert(phi < 2*twopi && phi > -twopi);
  double ptm = (m == 0) ? pt : sqrt(pt*pt+m*m);
  double exprap = exp(y);
  double pminus = ptm/exprap;
  double pplus  = ptm*exprap;
  double px = pt*cos(phi);
  double py = pt*sin(phi);
  PseudoJet mom(px,py,0.5*(pplus-pminus),0.5*(pplus+pminus));
  mom.set_cached_rap_phi(y,phi);
  return mom;
  //return PseudoJet(pt*cos(phi), pt*sin(phi), ptm*sinh(y), ptm*cosh(y));
}


//----------------------------------------------------------------------
// return kt-distance between this jet and another one
double PseudoJet::kt_distance(const PseudoJet & other) const {
  //double distance = min(this->kt2(), other.kt2());
  double distance = min(_kt2, other._kt2);
  double dphi = abs(phi() - other.phi());
  if (dphi > pi) {dphi = twopi - dphi;}
  double drap = rap() - other.rap();
  distance = distance * (dphi*dphi + drap*drap);
  return distance;
}


//----------------------------------------------------------------------
// return squared cylinder (eta-phi) distance between this jet and another one
double PseudoJet::plain_distance(const PseudoJet & other) const {
  double dphi = abs(phi() - other.phi());
  if (dphi > pi) {dphi = twopi - dphi;}
  double drap = rap() - other.rap();
  return (dphi*dphi + drap*drap);
}

//----------------------------------------------------------------------
/// returns other.phi() - this.phi(), i.e. the phi distance to
/// other, constrained to be in range -pi .. pi
double PseudoJet::delta_phi_to(const PseudoJet & other) const {
  double dphi = other.phi() - phi();
  if (dphi >  pi) dphi -= twopi;
  if (dphi < -pi) dphi += twopi;
  return dphi;
}


string PseudoJet::description() const{
  // the "default" case of a PJ which does not belong to any cluster sequence
  if (!_structure())
    return "standard PseudoJet (with no associated clustering information)";
  
  // for all the other cases, the description comes from the structure
  return _structure()->description();
}



//----------------------------------------------------------------------
//
// The following methods access the associated jet structure (if any)
//
//----------------------------------------------------------------------


//----------------------------------------------------------------------
// check whether this PseudoJet has an associated parent
// ClusterSequence
bool PseudoJet::has_associated_cluster_sequence() const{
  return (_structure()) && (_structure->has_associated_cluster_sequence());
}

//----------------------------------------------------------------------
// get a (const) pointer to the associated ClusterSequence (NULL if
// inexistent)
const ClusterSequence* PseudoJet::associated_cluster_sequence() const{
  if (! has_associated_cluster_sequence()) return NULL;

  return _structure->associated_cluster_sequence();
}


//----------------------------------------------------------------------
// check whether this PseudoJet has an associated parent
// ClusterSequence that is still valid
bool PseudoJet::has_valid_cluster_sequence() const{
  return (_structure()) && (_structure->has_valid_cluster_sequence());
}

//----------------------------------------------------------------------
// If there is a valid cluster sequence associated with this jet,
// returns a pointer to it; otherwise throws an Error.
//
// Open question: should these errors be upgraded to classes of their
// own so that they can be caught? [Maybe, but later]
const ClusterSequence * PseudoJet::validated_cs() const {
  return validated_structure_ptr()->validated_cs();
}


//----------------------------------------------------------------------
// set the associated structure
void PseudoJet::set_structure_shared_ptr(const SharedPtr<PseudoJetStructureBase> &structure_in){
  _structure = structure_in;
}

//----------------------------------------------------------------------
// return true if there is some strusture associated with this PseudoJet
bool PseudoJet::has_structure() const{
  return _structure();
}

//----------------------------------------------------------------------
// return a pointer to the structure (of type
// PseudoJetStructureBase*) associated with this PseudoJet.
//
// return NULL if there is no associated structure
const PseudoJetStructureBase* PseudoJet::structure_ptr() const {
  if (!_structure()) return NULL;
  return _structure();
}
  
//----------------------------------------------------------------------
// return a non-const pointer to the structure (of type
// PseudoJetStructureBase*) associated with this PseudoJet.
//
// return NULL if there is no associated structure
//
// Only use this if you know what you are doing. In any case,
// prefer the 'structure_ptr()' (the const version) to this method,
// unless you really need a write access to the PseudoJet's
// underlying structure.
PseudoJetStructureBase* PseudoJet::structure_non_const_ptr(){
  if (!_structure()) return NULL;
  return _structure();
}
  
//----------------------------------------------------------------------
// return a pointer to the structure (of type
// PseudoJetStructureBase*) associated with this PseudoJet.
//
// throw an error if there is no associated structure
const PseudoJetStructureBase* PseudoJet::validated_structure_ptr() const {
  if (!_structure()) 
    throw Error("Trying to access the structure of a PseudoJet which has no associated structure");
  return _structure();
}
  
//----------------------------------------------------------------------
// return a reference to the shared pointer to the
// PseudoJetStructureBase associated with this PseudoJet
const SharedPtr<PseudoJetStructureBase> & PseudoJet::structure_shared_ptr() const {
  return _structure;
}


//----------------------------------------------------------------------
// check if it has been recombined with another PseudoJet in which
// case, return its partner through the argument. Otherwise,
// 'partner' is set to 0.
//
// false is also returned if this PseudoJet has no associated
// ClusterSequence
bool PseudoJet::has_partner(PseudoJet &partner) const{
  return validated_structure_ptr()->has_partner(*this, partner);
}

//----------------------------------------------------------------------
// check if it has been recombined with another PseudoJet in which
// case, return its child through the argument. Otherwise, 'child'
// is set to 0.
// 
// false is also returned if this PseudoJet has no associated
// ClusterSequence, with the child set to 0
bool PseudoJet::has_child(PseudoJet &child) const{
  return validated_structure_ptr()->has_child(*this, child);
}

//----------------------------------------------------------------------
// check if it is the product of a recombination, in which case
// return the 2 parents through the 'parent1' and 'parent2'
// arguments. Otherwise, set these to 0.
//
// false is also returned if this PseudoJet has no parent
// ClusterSequence
bool PseudoJet::has_parents(PseudoJet &parent1, PseudoJet &parent2) const{
  return validated_structure_ptr()->has_parents(*this, parent1, parent2);
}

//----------------------------------------------------------------------
// check if the current PseudoJet contains the one passed as
// argument
//
// false is also returned if this PseudoJet has no associated
// ClusterSequence.
bool PseudoJet::contains(const PseudoJet &constituent) const{
  return validated_structure_ptr()->object_in_jet(constituent, *this);
}

//----------------------------------------------------------------------
// check if the current PseudoJet is contained the one passed as
// argument
//
// false is also returned if this PseudoJet has no associated
// ClusterSequence
bool PseudoJet::is_inside(const PseudoJet &jet) const{
  return validated_structure_ptr()->object_in_jet(*this, jet);
}


//----------------------------------------------------------------------
// returns true if the PseudoJet has constituents
bool PseudoJet::has_constituents() const{
  return (_structure()) && (_structure->has_constituents());
}

//----------------------------------------------------------------------
// retrieve the constituents.
vector<PseudoJet> PseudoJet::constituents() const{
  return validated_structure_ptr()->constituents(*this);
}


//----------------------------------------------------------------------
// returns true if the PseudoJet has support for exclusive subjets
bool PseudoJet::has_exclusive_subjets() const{
  return (_structure()) && (_structure->has_exclusive_subjets());
}

//----------------------------------------------------------------------
// return a vector of all subjets of the current jet (in the sense
// of the exclusive algorithm) that would be obtained when running
// the algorithm with the given dcut. 
//
// Time taken is O(m ln m), where m is the number of subjets that
// are found. If m gets to be of order of the total number of
// constituents in the jet, this could be substantially slower than
// just getting that list of constituents.
//
// an Error is thrown if this PseudoJet has no currently valid
// associated ClusterSequence
std::vector<PseudoJet> PseudoJet::exclusive_subjets (const double dcut) const {
  return validated_structure_ptr()->exclusive_subjets(*this, dcut);
}

//----------------------------------------------------------------------
// return the size of exclusive_subjets(...); still n ln n with same
// coefficient, but marginally more efficient than manually taking
// exclusive_subjets.size()
//
// an Error is thrown if this PseudoJet has no currently valid
// associated ClusterSequence
int PseudoJet::n_exclusive_subjets(const double dcut) const {
  return validated_structure_ptr()->n_exclusive_subjets(*this, dcut);
}

//----------------------------------------------------------------------
// return the list of subjets obtained by unclustering the supplied
// jet down to n subjets (or all constituents if there are fewer
// than n).
//
// requires n ln n time
//
// an Error is thrown if this PseudoJet has no currently valid
// associated ClusterSequence
std::vector<PseudoJet> PseudoJet::exclusive_subjets_up_to (int nsub) const {
  return validated_structure_ptr()->exclusive_subjets_up_to(*this, nsub);
}

//----------------------------------------------------------------------
// Same as exclusive_subjets_up_to but throws an error if there are
// fewer than nsub particles in the jet
std::vector<PseudoJet> PseudoJet::exclusive_subjets (int nsub) const {
  vector<PseudoJet> subjets = exclusive_subjets_up_to(nsub);
  if (int(subjets.size()) < nsub) {
    ostringstream err;
    err << "Requested " << nsub << " exclusive subjets, but there were only " 
	<< subjets.size() << " particles in the jet";
    throw Error(err.str());
  }
  return subjets;
}

//----------------------------------------------------------------------
// return the dij that was present in the merging nsub+1 -> nsub 
// subjets inside this jet.
//
// an Error is thrown if this PseudoJet has no currently valid
// associated ClusterSequence
double PseudoJet::exclusive_subdmerge(int nsub) const {
  return validated_structure_ptr()->exclusive_subdmerge(*this, nsub);
}

//----------------------------------------------------------------------
// return the maximum dij that occurred in the whole event at the
// stage that the nsub+1 -> nsub merge of subjets occurred inside 
// this jet.
//
// an Error is thrown if this PseudoJet has no currently valid
// associated ClusterSequence
double PseudoJet::exclusive_subdmerge_max(int nsub) const {
  return validated_structure_ptr()->exclusive_subdmerge_max(*this, nsub);
}


// returns true if a jet has pieces
//
// By default a single particle or a jet coming from a
// ClusterSequence have no pieces and this methos will return false.
bool PseudoJet::has_pieces() const{
  return ((_structure()) && (_structure->has_pieces(*this)));
}

// retrieve the pieces that make up the jet. 
//
// By default a jet does not have pieces.
// If the underlying interface supports "pieces" retrieve the
// pieces from there.
std::vector<PseudoJet> PseudoJet::pieces() const{
  return validated_structure_ptr()->pieces(*this);
  // if (!has_pieces())
  //   throw Error("Trying to retrieve the pieces of a PseudoJet that has no support for pieces.");
  //
  // return _structure->pieces(*this);
}


//----------------------------------------------------------------------
// the following ones require a computation of the area in the
// associated ClusterSequence (See ClusterSequenceAreaBase for details)
//----------------------------------------------------------------------

#ifndef __FJCORE__

//----------------------------------------------------------------------
// if possible, return a valid ClusterSequenceAreaBase pointer; otherwise
// throw an error
const ClusterSequenceAreaBase * PseudoJet::validated_csab() const {
  const ClusterSequenceAreaBase *csab = dynamic_cast<const ClusterSequenceAreaBase*>(validated_cs());
  if (csab == NULL) throw Error("you requested jet-area related information, but the PseudoJet does not have associated area information.");
  return csab;
}


//----------------------------------------------------------------------
// check if it has a defined area
bool PseudoJet::has_area() const{
  //if (! has_associated_cluster_sequence()) return false;
  if (! has_structure()) return false;
  return (validated_structure_ptr()->has_area() != 0);
}

//----------------------------------------------------------------------
// return the jet (scalar) area.
// throw an Error if there is no support for area in the associated CS
double PseudoJet::area() const{
  return validated_structure_ptr()->area(*this);
}

//----------------------------------------------------------------------
// return the error (uncertainty) associated with the determination
// of the area of this jet.
// throws an Error if there is no support for area in the associated CS
double PseudoJet::area_error() const{
  return validated_structure_ptr()->area_error(*this);
}

//----------------------------------------------------------------------
// return the jet 4-vector area
// throws an Error if there is no support for area in the associated CS
PseudoJet PseudoJet::area_4vector() const{
  return validated_structure_ptr()->area_4vector(*this);
}

//----------------------------------------------------------------------
// true if this jet is made exclusively of ghosts
// throws an Error if there is no support for area in the associated CS
bool PseudoJet::is_pure_ghost() const{
  return validated_structure_ptr()->is_pure_ghost(*this);
}

#endif  // __FJCORE__

//----------------------------------------------------------------------
//
// end of the methods accessing the information in the associated
// Cluster Sequence
//
//----------------------------------------------------------------------

//----------------------------------------------------------------------
/// provide a meaningful error message for InexistentUserInfo
PseudoJet::InexistentUserInfo::InexistentUserInfo() : Error("you attempted to perform a dynamic cast of a PseudoJet's extra info, but the extra info pointer was null")
{}


//----------------------------------------------------------------------
// sort the indices so that values[indices[0..n-1]] is sorted
// into increasing order 
void sort_indices(vector<int> & indices, 
			 const vector<double> & values) {
  IndexedSortHelper index_sort_helper(&values);
  sort(indices.begin(), indices.end(), index_sort_helper);
}



//----------------------------------------------------------------------
/// given a vector of values with a one-to-one correspondence with the
/// vector of objects, sort objects into an order such that the
/// associated values would be in increasing order
template<class T> vector<T>  objects_sorted_by_values(
                       const vector<T> & objects, 
		       const vector<double> & values) {

  assert(objects.size() == values.size());

  // get a vector of indices
  vector<int> indices(values.size());
  for (size_t i = 0; i < indices.size(); i++) {indices[i] = i;}
  
  // sort the indices
  sort_indices(indices, values);
  
  // copy the objects 
  vector<T> objects_sorted(objects.size());
  
  // place the objects in the correct order
  for (size_t i = 0; i < indices.size(); i++) {
    objects_sorted[i] = objects[indices[i]];
  }

  return objects_sorted;
}

//----------------------------------------------------------------------
/// return a vector of jets sorted into decreasing kt2
vector<PseudoJet> sorted_by_pt(const vector<PseudoJet> & jets) {
  vector<double> minus_kt2(jets.size());
  for (size_t i = 0; i < jets.size(); i++) {minus_kt2[i] = -jets[i].kt2();}
  return objects_sorted_by_values(jets, minus_kt2);
}

//----------------------------------------------------------------------
/// return a vector of jets sorted into increasing rapidity
vector<PseudoJet> sorted_by_rapidity(const vector<PseudoJet> & jets) {
  vector<double> rapidities(jets.size());
  for (size_t i = 0; i < jets.size(); i++) {rapidities[i] = jets[i].rap();}
  return objects_sorted_by_values(jets, rapidities);
}

//----------------------------------------------------------------------
/// return a vector of jets sorted into decreasing energy
vector<PseudoJet> sorted_by_E(const vector<PseudoJet> & jets) {
  vector<double> energies(jets.size());
  for (size_t i = 0; i < jets.size(); i++) {energies[i] = -jets[i].E();}
  return objects_sorted_by_values(jets, energies);
}

//----------------------------------------------------------------------
/// return a vector of jets sorted into increasing pz
vector<PseudoJet> sorted_by_pz(const vector<PseudoJet> & jets) {
  vector<double> pz(jets.size());
  for (size_t i = 0; i < jets.size(); i++) {pz[i] = jets[i].pz();}
  return objects_sorted_by_values(jets, pz);
}



//-------------------------------------------------------------------------------
// helper functions to build a jet made of pieces
//-------------------------------------------------------------------------------

// build a "CompositeJet" from the vector of its pieces
//
// In this case, E-scheme recombination is assumed to compute the
// total momentum
PseudoJet join(const vector<PseudoJet> & pieces){
  // compute the total momentum
  //--------------------------------------------------
  PseudoJet result;  // automatically initialised to 0
  for (unsigned int i=0; i<pieces.size(); i++)
    result += pieces[i];

  // attach a CompositeJetStructure to the result
  //--------------------------------------------------
  CompositeJetStructure *cj_struct = new CompositeJetStructure(pieces);

  result.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(cj_struct));

  return result;
}

// build a "CompositeJet" from a single PseudoJet
PseudoJet join(const PseudoJet & j1){
  return join(vector<PseudoJet>(1,j1));
}

// build a "CompositeJet" from two PseudoJet
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2){
  vector<PseudoJet> pieces;
  pieces.reserve(2);
  pieces.push_back(j1);
  pieces.push_back(j2);
  return join(pieces);
}

// build a "CompositeJet" from 3 PseudoJet
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3){
  vector<PseudoJet> pieces;
  pieces.reserve(3);
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  return join(pieces);
}

// build a "CompositeJet" from 4 PseudoJet
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3, const PseudoJet & j4){
  vector<PseudoJet> pieces;
  pieces.reserve(4);
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  pieces.push_back(j4);
  return join(pieces);
}




FASTJET_END_NAMESPACE

