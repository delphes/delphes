//FJSTARTHEADER
// $Id: RestFrameNSubjettinessTagger.cc 3433 2014-07-23 08:17:03Z salam $
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

#include <fastjet/tools/RestFrameNSubjettinessTagger.hh>
#include <fastjet/tools/Boost.hh>
#include <fastjet/ClusterSequence.hh>
#include <sstream>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//------------------------------------------------------------------------
// RestFrameNSubjettinessTagger class implementation
//------------------------------------------------------------------------

//------------------------------------------------------------------------
// tagger description
string RestFrameNSubjettinessTagger::description() const{ 
  ostringstream oss;
  oss << "RestFrameNSubjettiness tagger that performs clustering in the jet rest frame with " 
      << _subjet_def.description() 
      << ", supplemented with cuts tau_2 < " << _t2cut 
      << " and cos(theta_s) < " << _costscut;
  return oss.str();
}


//------------------------------------------------------------------------
// action on a single jet
// returns the tagged PseudoJet if successful, 0 otherwise
PseudoJet RestFrameNSubjettinessTagger::result(const PseudoJet & jet) const{
  // make sure that the jet has constituents
  if (!jet.has_constituents())
    throw("The jet you try to tag needs to have accessible constituents");
   
  // get the constituents and boost them into the rest frame of the jet
  vector<PseudoJet> rest_input = jet.constituents();
  for (unsigned int i=0; i<rest_input.size(); i++)
    rest_input[i].unboost(jet);

  ClusterSequence cs_rest(rest_input, _subjet_def);
  vector<PseudoJet> subjets = (_use_exclusive)
    ? cs_rest.exclusive_jets(2)
    : sorted_by_E(cs_rest.inclusive_jets());

  // impose the cuts in the rest-frame
  if (subjets.size()<2) return PseudoJet();

  const PseudoJet &j0 = subjets[0];
  const PseudoJet &j1 = subjets[1];

  /// impose the cut on cos(theta_s)
  double ct0 = (j0.px()*jet.px() + j0.py()*jet.py() + j0.pz()*jet.pz())
    /sqrt(j0.modp2()*jet.modp2());
  double ct1 = (j1.px()*jet.px() + j1.py()*jet.py() + j1.pz()*jet.pz())
    /sqrt(j1.modp2()*jet.modp2());
  if ((ct0 > _costscut) || (ct1 > _costscut)) return PseudoJet();
  
  // ccompute the 2-subjettiness and impose the coresponding cut
  double tau2 = 0.0;
  for (unsigned int i=0; i<rest_input.size(); i++)
    tau2 += min(dot_product(rest_input[i], j0), 
                dot_product(rest_input[i], j1));

  tau2 *= (2.0/jet.m2());

  if (tau2 > _t2cut) return PseudoJet();

  // We have a positive tag, 
  //  - boost everything back into the lab frame
  //  - record the info in the interface
  // Note that in order to point to the correct Clustersequence, the
  // subjets must be taken from the boosted one. We extract that
  // through the history index of the rest-frame subjets
  ClusterSequence * cs_structure = new ClusterSequence();
  Boost boost(jet);
  cs_structure->transfer_from_sequence(cs_rest, &boost);
  PseudoJet subjet_lab1 = cs_structure->jets()[cs_rest.history()[subjets[0].cluster_hist_index()].jetp_index];
  PseudoJet subjet_lab2 = cs_structure->jets()[cs_rest.history()[subjets[0].cluster_hist_index()].jetp_index];
    
  PseudoJet result_local = join<StructureType>(subjet_lab1,subjet_lab2);
  StructureType * s = (StructureType *) result_local.structure_non_const_ptr();
//  s->_original_jet = jet;
  s->_tau2 = tau2;
  s->_costhetas = max(ct0, ct1);

  // keep the rest-frame CS alive
  cs_structure->delete_self_when_unused();

  return result_local;
}


FASTJET_END_NAMESPACE
