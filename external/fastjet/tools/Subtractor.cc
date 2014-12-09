//FJSTARTHEADER
// $Id: Subtractor.cc 3670 2014-09-08 14:17:59Z soyez $
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

#include "fastjet/tools/Subtractor.hh"
#include <cassert>
#include <sstream>
#include <limits>
using namespace std;

FASTJET_BEGIN_NAMESPACE     // defined in fastjet/internal/base.hh

const double Subtractor::_invalid_rho = -numeric_limits<double>::infinity();


//----------------------------------------------------------------------
// ctor
Subtractor::Subtractor(double rho) : _bge(0), _rho(rho) {
  if (_rho<0.0) throw Error("Subtractor(rho) was passed a negative rho value; rho should be >= 0");
  set_defaults();
}

//----------------------------------------------------------------------
// ctor
Subtractor::Subtractor(double rho, double rho_m) : _bge(0), _rho(rho) {
  if (_rho<0.0) throw Error("Subtractor(rho, rho_m) was passed a negative rho value; rho should be >= 0");
  if (rho_m<0.0) throw Error("Subtractor(rho, rho_m) was passed a negative rho_m value; rho_m should be >= 0");
  set_defaults();
  _rho_m = rho_m;
  set_use_rho_m(true);
}

//----------------------------------------------------------------------
void Subtractor::set_defaults(){
  _rho_m = _invalid_rho;
  _use_rho_m = false; // likely to change in future releases!!
  _safe_mass = false; // likely to change in future releases!!

  _sel_known_vertex = Selector();
  _sel_leading_vertex = Selector();
}

//----------------------------------------------------------------------
// perform the subtraction of a given jet
PseudoJet Subtractor::result(const PseudoJet & jet) const {
  if (!jet.has_area()){
    throw Error("Subtractor::result(...): Trying to subtract a jet without area support");
  }

  PseudoJet known_lv, known_pu;
  PseudoJet unknown = jet;
  if (_sel_known_vertex.worker()){
    // separate the jet constituents in 3 groups:
    //   unknown vertex
    //   known vertex, leading vertex
    //   known vertex, non-leading vertex (PU)
    vector<PseudoJet> constits_unknown, constits_known;
    _sel_known_vertex.sift(jet.constituents(), 
			   constits_known,
			   constits_unknown);
    vector<PseudoJet> constits_known_lv, constits_known_pu;
    _sel_leading_vertex.sift(constits_known,
			     constits_known_lv,
			     constits_known_pu);

    // For the parts related to the known vertices (LV or PU), we just
    // sum the 4-momenta. For the unknown part, we assign it the full
    // jet area.
    known_lv = (constits_known_lv.size()!=0) 
      ? SelectorIdentity().sum(constits_known_lv) : 0.0*jet;
    known_pu = (constits_known_pu.size()!=0) 
      ? SelectorIdentity().sum(constits_known_pu) : 0.0*jet;
    if (constits_unknown.size()==0){
      // no need for any form of subtraction!
      PseudoJet subtracted_jet = jet;
      subtracted_jet.reset_momentum(known_lv);
      return subtracted_jet;
    }
    unknown = jet; // that keeps all info including area
    unknown.reset_momentum(SelectorIdentity().sum(constits_unknown));
  } else {
    known_lv = jet; // ensures correct rap-phi!
    known_lv *= 0.0;
    known_pu = known_lv;
  }

  // prepare for the subtraction and compute the 4-vector to be
  // subtracted
  PseudoJet subtracted_jet = jet;
  PseudoJet to_subtract = known_pu + _amount_to_subtract(unknown);

  // sanity check for the transverse momentum
  if (to_subtract.pt2() < jet.pt2() ) { 
    // this subtraction should retain the jet's structural
    // information
    subtracted_jet -= to_subtract;
  } else { 
    // this sets the jet's momentum while maintaining all of the jet's
    // structural information
    subtracted_jet.reset_momentum(known_lv);
    return subtracted_jet;
  }

  // make sure that in the end the pt is at least the one known to
  // come from the leading vertex
  if (subtracted_jet.pt2() < known_lv.pt2()){
    subtracted_jet.reset_momentum(known_lv);
    return subtracted_jet;
  }

  // sanity check for the mass (if needed)
  if ((_safe_mass) && (subtracted_jet.m2() < known_lv.m2())){
    // in this case, we keep pt and phi as obtained from the
    // subtraction above and take rap and m from the part that comes
    // from the leading vertex (or the original jet if nothing comes
    // from the leading vertex)
    subtracted_jet.reset_momentum(PtYPhiM(subtracted_jet.pt(),
					  known_lv.rap(),
					  subtracted_jet.phi(),
					  known_lv.m()));
  }

  return subtracted_jet;
}

//----------------------------------------------------------------------
std::string Subtractor::description() const{
  if (_bge != 0) {
    string desc = "Subtractor that uses the following background estimator to determine rho: "+_bge->description();
    if (use_rho_m()) desc += "; including the rho_m correction";
    if (safe_mass()) desc += "; including mass safety tests";
    if (_sel_known_vertex.worker()){
      desc += "; using known vertex selection: "+_sel_known_vertex.description()+" and leading vertex selection: "+_sel_leading_vertex.description();
    }
    return desc;
  } else if (_rho != _invalid_rho) {
    ostringstream ostr;
    ostr << "Subtractor that uses a fixed value of rho = " << _rho;
    if (use_rho_m()) ostr << " and rho_m = " << _rho_m;
    return ostr.str();
  } else {
    return "Uninitialised subtractor";
  }
}

//----------------------------------------------------------------------
// compute the 4-vector that should be subtracted from the given
// jet
PseudoJet Subtractor::_amount_to_subtract(const PseudoJet &jet) const{
  // the "transverse momentum" part
  double rho;
  if (_bge != 0) {
    rho = _bge->rho(jet);
  } else if (_rho != _invalid_rho) {
    rho = _rho;
  } else {
    throw Error("Subtractor::_amount_to_subtract(...): default Subtractor does not have any information about the background, needed to perform the subtraction");
  }

  PseudoJet area = jet.area_4vector();
  PseudoJet to_subtract = rho*area;

  double const rho_m_warning_threshold = 1e-5;

  // add an optional contribution from the unknown particles masses
  if (_use_rho_m) {
    double rho_m;
    
    if (_bge != 0) {
      if (!_bge->has_rho_m()) throw Error("Subtractor::_amount_to_subtract(...): requested subtraction with rho_m from a background estimator, but the estimator does not have rho_m support");
      rho_m = _bge->rho_m(jet);
    } else if (_rho_m != _invalid_rho) {
      rho_m = _rho_m;
    } else {
      throw Error("Subtractor::_amount_to_subtract(...): default Subtractor does not have any information about the background rho_m, needed to perform the rho_m subtraction");
    }
    to_subtract += rho_m * PseudoJet(0.0, 0.0, area.pz(), area.E());
  } else if (_bge && 
             _bge->has_rho_m() && 
             _bge->rho_m(jet) > rho_m_warning_threshold * rho) {
    _unused_rho_m_warning.warn("Subtractor::_amount_to_subtract(...): Background estimator indicates non-zero rho_m, but use_rho_m()==false in subtractor; consider calling set_use_rho_m(true) to include the rho_m information");
  }

  return to_subtract;
}


FASTJET_END_NAMESPACE
