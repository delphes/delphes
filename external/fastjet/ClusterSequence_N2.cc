//FJSTARTHEADER
// $Id$
//
// Copyright (c) 2005-2025, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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


// The plain N^2 part of the ClusterSequence class -- separated out
// from the rest of the class implementation so as to speed up
// compilation of this particular part while it is under test.

#include "fastjet/internal/ClusterSequence_N2.icc"

#include<iostream>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


using namespace std;


//*************************************************************************
//
//                             THINGS FOR E+E-
//
//*************************************************************************


//----------------------------------------------------------------------
template<> inline void ClusterSequence::_bj_set_jetinfo(
                           EEBriefJet * const jetA, const int _jets_index) const {

  double E = _jets[_jets_index].E();
  double scale = E*E; // the default energy scale for the kt alg
  double p  = jet_def().extra_param(); // in case we're ee_genkt
  switch (_jet_algorithm) {
  case ee_kt_algorithm:
    assert(_Rparam > 2.0); // force this to be true! [not best place, but works]
    // recall that _invR2 is artificially set to 1 for this alg
    // so that we automatically have dij = scale * 2(1-cos theta_ij)
    // Normally, _Rparam should be automatically set to 4 from JetDefinition
    break; 
  case ee_genkt_algorithm:
    if (p <= 0 && scale < 1e-300) scale = 1e-300; // same dodgy safety as genkt
    scale = pow(scale,p);
    break;
  default:
    throw Error("Unrecognised jet algorithm");
  }
  jetA->kt2  = scale; // "kt2" might one day be renamed as "scale" or some such

  double norm = _jets[_jets_index].modp2();
  if (norm > 0) {
    norm = 1.0/sqrt(norm);
    jetA->nx = norm * _jets[_jets_index].px();
    jetA->ny = norm * _jets[_jets_index].py();
    jetA->nz = norm * _jets[_jets_index].pz();
  } else {
    jetA->nx = 0.0;
    jetA->ny = 0.0;
    jetA->nz = 1.0;
  }
  jetA->_jets_index = _jets_index;
  // initialise NN info as well
  jetA->NN_dist = _R2;
  jetA->NN      = NULL;
}

// this is declared to ensure that calls with the EEAccurateBriefJet
// are redirected to the EEBriefJet implementation (otherwise
// it tries the default pp template)
template<> inline void ClusterSequence::_bj_set_jetinfo(
                           EEAccurateBriefJet * const jetA, const int _jets_index) const {
  _bj_set_jetinfo<EEBriefJet>(jetA, _jets_index);
}


//----------------------------------------------------------------------
// returns the angular distance between the two jets, defined as
// 2*(1-cos theta_ab)
template<> double ClusterSequence::_bj_dist(
                const EEBriefJet * const jeta, 
                const EEBriefJet * const jetb) const {
  double dist = 1.0 
    - jeta->nx*jetb->nx
    - jeta->ny*jetb->ny
    - jeta->nz*jetb->nz;

  return dist*2; // distance is _2_*min(Ei^2,Ej^2)*(1-cos theta)
}

//----------------------------------------------------------------------
// returns the angular distance between the two jets, defined as
// 2*(1-cos theta_ab)
template<> double ClusterSequence::_bj_dist(
                const EEAccurateBriefJet * const jeta, 
                const EEAccurateBriefJet * const jetb) const {
  double dist = 1.0 
    - jeta->nx*jetb->nx
    - jeta->ny*jetb->ny
    - jeta->nz*jetb->nz;

  // if the distance is smaller than sqrt(epsilon), then switch to
  // cross-product based evaluation; this is intended to ensure an
  // absolute accuracy on sqrt(bj_dist), that is somewhere in the region
  // of max(epsilon, sqrt(epsilon) * dist) rather than epsilon.
  //
  // Let's write cos(theta) == cos and sin(theta)==sin; then we have
  //
  //    2*(1-cos) = 2*(1-cos^2)/(1+cos) = 2*sin^2/(1+cos);   
  //
  // To save a division, we then replace 2/(1+cos) -> 1, which introduces
  // a relative error of order sin^2(theta), which is sqrt(epsilon) when
  // dist is itself sqrt(epsilon). For that value of dist, the relative
  // error on the normal dot-product calculation of dist is itself
  // sqrt(epsilon). This motivates the choice switchover point.
  //
  // This approach has been adopted from PanScales work.
  if (dist*dist < numeric_limits<double>::epsilon()) {
    double cross_x = jeta->ny * jetb->nz - jetb->ny * jeta->nz;
    double cross_y = jeta->nz * jetb->nx - jetb->nz * jeta->nx;
    double cross_z = jeta->nx * jetb->ny - jetb->nx * jeta->ny;
    
    // 2(1-cos(theta)) ~ theta^2, which is |cross_product|^2
    dist = cross_x*cross_x + cross_y*cross_y + cross_z*cross_z;
    return dist;
  }

  // this alternative code has higher accuracy in some boundary regions
  // but is almost a factor of two slower than the plain dot product
  // (the code above is only 15% slower on a standard-looking event).  
  //if (dist < 1) {
  //  double cross_x = jeta->ny * jetb->nz - jetb->ny * jeta->nz;
  //  double cross_y = jeta->nz * jetb->nx - jetb->nz * jeta->nx;
  //  double cross_z = jeta->nx * jetb->ny - jetb->nx * jeta->ny;
  //  
  //  // 2(1-cos(theta)) ~ theta^2, which is |cross_product|^2
  //  double sinsqr = cross_x*cross_x + cross_y*cross_y + cross_z*cross_z;
  //  return 2*sinsqr/(2.0-dist);
  //  //return dist;
  //}

  return dist*2; // distance is _2_*min(Ei^2,Ej^2)*(1-cos theta)
}



// get explicit copies of the two N2 cluster cases we need
// plain BriefJet
void ClusterSequence::_simple_N2_cluster_BriefJet() {  
  _simple_N2_cluster<BriefJet>();
}


// e+e- BriefJet
void ClusterSequence::_simple_N2_cluster_EEBriefJet() {  
  _simple_N2_cluster<EEBriefJet>();
}

void ClusterSequence::_simple_N2_cluster_EEAccurateBriefJet() {  
  _simple_N2_cluster<EEAccurateBriefJet>();
}


// //----------------------------------------------------------------------
// /// Force instantiation of desired versions of _simple_N2_cluster
// ///
// /// This is not very elegant...
// void ClusterSequence::_dummy_N2_cluster_instantiation() {
//   _simple_N2_cluster<BriefJet>();
//   _simple_N2_cluster<EEBriefJet>();
// }

FASTJET_END_NAMESPACE

