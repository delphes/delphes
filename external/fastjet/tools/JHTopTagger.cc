//FJSTARTHEADER
// $Id: JHTopTagger.cc 3433 2014-07-23 08:17:03Z salam $
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

#include <fastjet/tools/JHTopTagger.hh>
#include <fastjet/Error.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <sstream>
#include <limits>

FASTJET_BEGIN_NAMESPACE

using namespace std;

//----------------------------------------------------------------------
// JHTopTagger class implementation
//----------------------------------------------------------------------

LimitedWarning JHTopTagger::_warnings_nonca;

//------------------------------------------------------------------------
// description of the tagger
string JHTopTagger::description() const{ 
  ostringstream oss;
  oss << "JHTopTagger with delta_p=" << _delta_p << ", delta_r=" << _delta_r
      << ", cos_theta_W_max=" << _cos_theta_W_max
      << " and mW = " << _mW;
  oss << description_of_selectors();
  return oss.str();
}

//------------------------------------------------------------------------
// returns the tagged PseudoJet if successful, 0 otherwise
//  - jet   the PseudoJet to tag
PseudoJet JHTopTagger::result(const PseudoJet & jet) const{
  // make sure that there is a "regular" cluster sequence associated
  // with the jet. Note that we also check it is valid (to avoid a
  // more criptic error later on)
  if (!jet.has_valid_cluster_sequence()){
    throw Error("JHTopTagger can only be applied on jets having an associated (and valid) ClusterSequence");
  }

  // warn if the jet has not been clustered with a Cambridge/Aachen
  // algorithm
  if (jet.validated_cs()->jet_def().jet_algorithm() != cambridge_algorithm)
    _warnings_nonca.warn("JHTopTagger should only be applied on jets from a Cambridge/Aachen clustering; use it with other algorithms at your own risk.");


  // do the first splitting
  vector<PseudoJet> split0 = _split_once(jet, jet);
  if (split0.size() == 0) return PseudoJet();

  // now try a second splitting on each of the resulting objects
  vector<PseudoJet> subjets;
  for (unsigned i = 0; i < 2; i++) {
    vector<PseudoJet> split1 = _split_once(split0[i], jet);
    if (split1.size() > 0) {
      subjets.push_back(split1[0]);
      subjets.push_back(split1[1]);
    } else {
      subjets.push_back(split0[i]);
    }
  }

  // make sure things make sense
  if (subjets.size() < 3) return PseudoJet();

  // now find the pair of objects closest in mass to the W
  double dmW_min = numeric_limits<double>::max();
  int ii=-1, jj=-1;
  for (unsigned i = 0 ; i < subjets.size()-1; i++) {
    for (unsigned j = i+1 ; j < subjets.size(); j++) {
      double dmW = abs(_mW - (subjets[i]+subjets[j]).m());
      if (dmW < dmW_min) {
        dmW_min = dmW; ii = i; jj = j;
      }
    }
  }

  // order the subjets in the following order:
  //  - hardest of the W subjets
  //  - softest of the W subjets
  //  - hardest of the remaining subjets
  //  - softest of the remaining subjets (if any)
  if (ii>0) std::swap(subjets[ii], subjets[0]);
  if (jj>1) std::swap(subjets[jj], subjets[1]);
  if (subjets[0].perp2() < subjets[1].perp2()) std::swap(subjets[0], subjets[1]);
  if ((subjets.size()>3) && (subjets[2].perp2() < subjets[3].perp2())) 
    std::swap(subjets[2], subjets[3]);
  
  // create the result and its structure
  const JetDefinition::Recombiner *rec
    = jet.associated_cluster_sequence()->jet_def().recombiner();

  PseudoJet W = join(subjets[0], subjets[1], *rec);
  PseudoJet non_W;
  if (subjets.size()>3) {
    non_W = join(subjets[2], subjets[3], *rec);
  } else {
    non_W = join(subjets[2], *rec);
  }
  PseudoJet result_local = join<JHTopTaggerStructure>(W, non_W, *rec);
  JHTopTaggerStructure *s = (JHTopTaggerStructure*) result_local.structure_non_const_ptr();
  s->_cos_theta_w = _cos_theta_W(result_local);

  // if the polarisation angle does not pass the cut, consider that
  // the tagging has failed
  //
  // Note that we could perhaps ensure this cut before constructing
  // the result structure but this has the advantage that the top
  // 4-vector is already available and does not have to de re-computed
  if (s->_cos_theta_w >= _cos_theta_W_max ||
      ! _top_selector.pass(result_local) || ! _W_selector.pass(W)
      ) {
    result_local *= 0.0;
  }

  return result_local;

  // // old version
  // PseudoJet result = join<JHTopTaggerStructure>(subjets, *rec);
  // JHTopTaggerStructure *s = (JHTopTaggerStructure*) result.structure_non_const_ptr();
  // //  s->_original_jet = jet;
  // s->_W = join(subjets[0], subjets[1], *rec);
  // if (subjets.size()>3)
  //   s->_non_W = join(subjets[2], subjets[3], *rec);
  // else
  //   s->_non_W = join(subjets[2], *rec);
  // s->_cos_theta_w = _cos_theta_W(result);
  // 
  // // if the polarisation angle does not pass the cut, consider that
  // // the tagging has failed
  // //
  // // Note that we could perhaps ensure this cut before constructing
  // // the result structure but this has the advantage that the top
  // // 4-vector is already available and does not have to de re-computed
  // if (s->_cos_theta_w >= _cos_theta_W_max)
  //   return PseudoJet();
  // 
  // return result;
}

// runs the Johns Hopkins decomposition procedure
vector<PseudoJet> JHTopTagger::_split_once(const PseudoJet & jet_to_split,
                                           const PseudoJet & reference_jet) const{
  PseudoJet this_jet = jet_to_split;
  PseudoJet p1, p2;
  vector<PseudoJet> result_local;
  while (this_jet.has_parents(p1, p2)) {
    if (p2.perp2() > p1.perp2()) std::swap(p1,p2); // order with hardness
    if (p1.perp() < _delta_p * reference_jet.perp()) break; // harder is too soft wrt original jet
    if ( (abs(p2.rap()-p1.rap()) + abs(p2.delta_phi_to(p1))) < _delta_r) break; // distance is too small
    if (p2.perp() < _delta_p * reference_jet.perp()) {
      this_jet = p1; // softer is too soft wrt original, so ignore it
      continue; 
    }
    //result.push_back(this_jet);
    result_local.push_back(p1);
    result_local.push_back(p2);
    break;
  }
  return result_local;
}




FASTJET_END_NAMESPACE

