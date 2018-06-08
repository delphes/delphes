// $Id: BottomUpSoftDrop.cc 1064 2017-09-08 09:19:57Z gsoyez $
//
// Copyright (c) 2017-, Gavin P. Salam, Gregory Soyez, Jesse Thaler,
// Kevin Zhou, Frederic Dreyer
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#include "BottomUpSoftDrop.hh"
#include <cassert>
#include <algorithm>
#include <sstream>
#include <typeinfo>
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/Selector.hh"
#include "fastjet/config.h"


using namespace std;


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//----------------------------------------------------------------------
// BottomUpSoftDrop class
//----------------------------------------------------------------------
  
// action on a single jet
PseudoJet BottomUpSoftDrop::result(const PseudoJet &jet) const{
  // soft drop can only be applied to jets that have constituents
  if (!jet.has_constituents()){
    throw Error("BottomUpSoftDrop: trying to apply the Soft Drop transformer to a jet that has no constituents");
  }

  // if the jet has area support and there are explicit ghosts, we can
  // transfer that support to the internal re-clustering
  //
  // Note that this is just meant to maintain the information since
  // all the jes will have a 0 area
  bool do_areas = jet.has_area() && _check_explicit_ghosts(jet);

  // build the soft drop plugin
  BottomUpSoftDropPlugin * softdrop_plugin;

  // for some constructors, we get the recombiner from the 
  // input jet -- some acrobatics are needed
  if (_get_recombiner_from_jet) {
    JetDefinition jet_def = _jet_def;

    // if all the pieces have a shared recombiner, we'll use that
    // one. Otherwise, use the one from _jet_def as a fallback.
    JetDefinition jet_def_for_recombiner;
    if (_check_common_recombiner(jet, jet_def_for_recombiner)){
#if FASTJET_VERSION_NUMBER >= 30100
      // Note that this is better than the option directly passing the
      // recombiner (for cases where th ejet def own its recombiner)
      // but it's only available in FJ>=3.1
      jet_def.set_recombiner(jet_def_for_recombiner);
#else
      jet_def.set_recombiner(jet_def_for_recombiner.recombiner());
#endif
    }
    softdrop_plugin = new BottomUpSoftDropPlugin(jet_def, _beta, _symmetry_cut, _R0);
  } else {
    softdrop_plugin = new BottomUpSoftDropPlugin(_jet_def, _beta, _symmetry_cut, _R0);
  }

  // now recluster the constituents of the jet with that plugin
  JetDefinition internal_jet_def(softdrop_plugin);
  // flag the plugin for automatic deletion _before_ we make
  // copies (so that as long as the copies are also present
  // it doesn't get deleted).
  internal_jet_def.delete_plugin_when_unused();

  ClusterSequence * cs;
  if (do_areas){
    vector<PseudoJet> particles, ghosts;
    SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
    // determine the ghost area from the 1st ghost (if none, any value
    // will do, as the area will be 0 and subtraction will have
    // no effect!)
    double ghost_area = (ghosts.size()) ? ghosts[0].area() : 0.01;
    cs = new ClusterSequenceActiveAreaExplicitGhosts(particles, internal_jet_def, 
  						     ghosts, ghost_area);
  } else {
    cs = new ClusterSequence(jet.constituents(), internal_jet_def);
  }
  
  PseudoJet result_local = SelectorNHardest(1)(cs->inclusive_jets())[0];
  BottomUpSoftDropStructure * s = new BottomUpSoftDropStructure(result_local);
  s->_beta = _beta;
  s->_symmetry_cut = _symmetry_cut;
  s->_R0   = _R0;
  result_local.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(s));
  
  // make sure things remain persistent -- i.e. tell the jet definition
  // and the cluster sequence that it is their responsibility to clean 
  // up memory once the "result" reaches the end of its life in the user's
  // code. (The CS deletes itself when the result goes out of scope and
  // that also triggers deletion of the plugin)
  cs->delete_self_when_unused();

  return result_local;  
}

// global grooming on a full event
// note: does not support jet areas
vector<PseudoJet> BottomUpSoftDrop::global_grooming(const vector<PseudoJet> & event) const {
  // start by reclustering the event into one very large jet
  ClusterSequence cs(event, _jet_def);
  std::vector<PseudoJet> global_jet = SelectorNHardest(1)(cs.inclusive_jets());
  // if the event is empty, do nothing
  if (global_jet.size() == 0) return vector<PseudoJet>();
  
  // apply the groomer to the large jet
  PseudoJet result = this->result(global_jet[0]);
  return result.constituents();
}
  
// check if the jet has explicit_ghosts (knowing that there is an
// area support)
bool BottomUpSoftDrop::_check_explicit_ghosts(const PseudoJet &jet) const {
  // if the jet comes from a Clustering check explicit ghosts in that
  // clustering
  if (jet.has_associated_cluster_sequence())
    return jet.validated_csab()->has_explicit_ghosts();

  // if the jet has pieces, recurse in the pieces
  if (jet.has_pieces()){
    vector<PseudoJet> pieces = jet.pieces();
    for (unsigned int i=0;i<pieces.size(); i++)
      if (!_check_explicit_ghosts(pieces[i])) return false;
    // never returned false, so we're OK.
    return true;
  }

  // return false for any other (unknown) structure
  return false;
}

// see if there is a common recombiner among the pieces; if there
// is return true and set jet_def_for_recombiner so that the
// recombiner can be taken from that JetDefinition. Otherwise,
// return false. 'assigned' is initially false; when true, each
// time we meet a new jet definition, we'll check it shares the
// same recombiner as jet_def_for_recombiner.
bool BottomUpSoftDrop::_check_common_recombiner(const PseudoJet &jet, 
						JetDefinition &jet_def_for_recombiner,
						bool assigned) const {
  if (jet.has_associated_cluster_sequence()){
    // if the jet def for recombination has already been assigned, check if we have the same
    if (assigned)
      return jet.validated_cs()->jet_def().has_same_recombiner(jet_def_for_recombiner);

    // otherwise, assign it.
    jet_def_for_recombiner = jet.validated_cs()->jet_def();
    assigned = true;
    return true;
  }

  // if the jet has pieces, recurse in the pieces
  if (jet.has_pieces()){
    vector<PseudoJet> pieces = jet.pieces();
    if (pieces.size() == 0) return false;
    for (unsigned int i=0;i<pieces.size(); i++)
      if (!_check_common_recombiner(pieces[i], jet_def_for_recombiner, assigned)) return false;
    // never returned false, so we're OK.
    return true;
  }
  // return false for any other (unknown) structure
  return false;
}


// description
string BottomUpSoftDrop::description() const{
  ostringstream oss;
  oss << "BottomUpSoftDrop with jet_definition = (" << _jet_def.description() << ")"
      << ", symmetry_cut = " << _symmetry_cut	<< ", beta = " << _beta << ", R0 = " << _R0;
  return oss.str();
}
  
//----------------------------------------------------------------------
// BottomUpSoftDropRecombiner class
//----------------------------------------------------------------------

// perform a recombination taking into account the Soft Drop
// conditions
void BottomUpSoftDropRecombiner::recombine(const PseudoJet &pa, 
                                           const PseudoJet &pb,
                                           PseudoJet &pab) const {
  PseudoJet p;
  _recombiner->recombine(pa, pb, p);

  // calculate symmetry parameter
  double symmetry_cut_fn =_symmetry_cut * pow(pa.squared_distance(pb)/_R0sqr, 0.5*_beta);
  double pta = pa.pt();
  double ptb = pb.pt();
  // make sure denominator is non-zero
  double sym = pta + ptb;
  if (sym == 0){
    pab = p;
    return;
  }
  sym = min(pta, ptb) / sym;
    
  // if the 2 particles pass the soft drop condition, do the
  // recombination
  if (sym > symmetry_cut_fn){
    pab=p;
    return;
  }

  // if the soft drop condition is not passed, return the particle
  // with largest pt
  if (pta < ptb) {
    pab = pb;
    _rejected.push_back(pa.cluster_hist_index());
  } else {
    pab = pa;
    _rejected.push_back(pb.cluster_hist_index());
  }
}

//----------------------------------------------------------------------
// BottomUpSoftDropPlugin class
//----------------------------------------------------------------------

// the actual clustering work for the plugin
void BottomUpSoftDropPlugin::run_clustering(ClusterSequence &input_cs) const {
  // declare a soft drop recombiner
  BottomUpSoftDropRecombiner softdrop_recombiner(_beta, _symmetry_cut, _R0, _jet_def.recombiner());
  JetDefinition jet_def = _jet_def;
  jet_def.set_recombiner(&softdrop_recombiner);

  // cluster the particles using that recombiner
  ClusterSequence internal_cs(input_cs.jets(), jet_def);
  const vector<ClusterSequence::history_element> & internal_hist = internal_cs.history();

  // transfer the list of "childless" elements into a bool vector
  vector<bool> kept(internal_hist.size(), true);
  const vector<unsigned int> &sd_rej = softdrop_recombiner.rejected();
  for (unsigned int i=0;i<sd_rej.size(); i++) kept[sd_rej[i]]=false;

  // browse the history, keeping only the elements that have not been
  // vetoed.
  //
  // In the process we build a map for the history indices
  vector<unsigned int> internal2input(internal_hist.size());
  for (unsigned int i=0; i<input_cs.jets().size(); i++)
    internal2input[i] = i;

  for (unsigned int i=input_cs.jets().size(); i<internal_hist.size(); i++){
    const ClusterSequence::history_element &he = internal_hist[i];

    // deal with recombinations with the beam
    if (he.parent2 == ClusterSequence::BeamJet){
      int internal_jetp_index = internal_hist[he.parent1].jetp_index;
      int internal_hist_index = internal_cs.jets()[internal_jetp_index].cluster_hist_index();

      int input_jetp_index = input_cs.history()[internal2input[internal_hist_index]].jetp_index;

      // cout << "Beam recomb for internal " << internal_hist_index
      //            << " (input jet index=" << input_jetp_index << endl;

      input_cs.plugin_record_iB_recombination(input_jetp_index, he.dij);
      continue;
    }

    // now, deal with two-body recombinations
    if (!kept[he.parent1]){ // 1 is rejected, we keep only 2
      internal2input[i]=internal2input[he.parent2];
      // cout << "rejecting internal " << he.parent1
      //            << ", mapping internal " << i 
      //            << " to internal " << he.parent2
      //            << " i.e. " << internal2input[i] << endl;
    } else if (!kept[he.parent2]){ // 2 is rejected, we keep only 1
      internal2input[i]=internal2input[he.parent1];
      // cout << "rejecting internal " << he.parent2 
      //            << ", mapping internal " << i 
      //            << " to internal " << he.parent1
      //            << " i.e. " << internal2input[i] << endl;
    } else { // do the recombination
      int new_index;
      input_cs.plugin_record_ij_recombination(input_cs.history()[internal2input[he.parent1]].jetp_index,
					      input_cs.history()[internal2input[he.parent2]].jetp_index,
					      he.dij, internal_cs.jets()[he.jetp_index], new_index);
      internal2input[i]=input_cs.jets()[new_index].cluster_hist_index();
      // cout << "merging " << internal2input[he.parent1] << " (int: " << he.parent1 << ")"
      //      << " and "    << internal2input[he.parent2] << " (int: " << he.parent2 << ")"
      //      << " into "   << internal2input[i] << " (int: " << i << ")" << endl;
    }
  }
}

// description of the plugin
string BottomUpSoftDropPlugin::description() const {
  ostringstream oss;
  oss << "BottomUpSoftDropPlugin with jet_definition = (" << _jet_def.description()
      <<"), symmetry_cut = " << _symmetry_cut << ", beta = "
      << _beta << ", R0 = " << _R0;
  return oss.str();
}


}

FASTJET_END_NAMESPACE
