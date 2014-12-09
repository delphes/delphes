//FJSTARTHEADER
// $Id: Pruner.cc 3481 2014-07-29 17:24:12Z soyez $
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

#include "fastjet/tools/Pruner.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/Selector.hh"
#include <cassert>
#include <algorithm>
#include <sstream>
#include <typeinfo>

using namespace std;


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


//----------------------------------------------------------------------
// class Pruner
//----------------------------------------------------------------------

//----------------------------------------------------------------------
// alternative (dynamic) ctor
//  \param jet_def     the jet definition for the internal clustering
//  \param zcut_dyn    dynamic pt-fraction cut in the pruning
//  \param Rcut_dyn    dynamic angular distance cut in the pruning
Pruner::Pruner(const JetDefinition &jet_def, 
         const FunctionOfPseudoJet<double> *zcut_dyn,
         const FunctionOfPseudoJet<double> *Rcut_dyn)
  : _jet_def(jet_def), _zcut(0), _Rcut_factor(0),
    _zcut_dyn(zcut_dyn), _Rcut_dyn(Rcut_dyn), _get_recombiner_from_jet(false)  {
  assert(_zcut_dyn != 0 && _Rcut_dyn != 0);
}

//----------------------------------------------------------------------
// action on a single jet
PseudoJet Pruner::result(const PseudoJet &jet) const{
  // pruning can only be applied to jets that have constituents
  if (!jet.has_constituents()){
    throw Error("Pruner: trying to apply the Pruner transformer to a jet that has no constituents");
  }

  // if the jet has area support and there are explicit ghosts, we can
  // transfer that support to the internal re-clustering
  bool do_areas = jet.has_area() && _check_explicit_ghosts(jet);

  // build the pruning plugin
  double Rcut = (_Rcut_dyn) ? (*_Rcut_dyn)(jet) : _Rcut_factor * 2.0*jet.m()/jet.perp();
  double zcut = (_zcut_dyn) ? (*_zcut_dyn)(jet) : _zcut;
  PruningPlugin * pruning_plugin;

  // for some constructors, we get the recombiner from the 
  // input jet -- some acrobatics are needed
  if (_get_recombiner_from_jet) {
    JetDefinition jet_def = _jet_def;

    // if all the pieces have a shared recombiner, we'll use that
    // one. Otherwise, use the one from _jet_def as a fallback.
    JetDefinition jet_def_for_recombiner;
    if (_check_common_recombiner(jet, jet_def_for_recombiner)){
      jet_def.set_recombiner(jet_def_for_recombiner);
    }
    pruning_plugin = new PruningPlugin(jet_def, zcut, Rcut);
  } else {
    pruning_plugin = new PruningPlugin(_jet_def, zcut, Rcut);
  }

  // now recluster the constituents of the jet with that plugin
  JetDefinition internal_jet_def(pruning_plugin);
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
  PrunerStructure * s = new PrunerStructure(result_local);
  s->_Rcut = Rcut;
  s->_zcut = zcut;
  result_local.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(s));
  
  // make sure things remain persistent -- i.e. tell the jet definition
  // and the cluster sequence that it is their responsibility to clean 
  // up memory once the "result" reaches the end of its life in the user's
  // code. (The CS deletes itself when the result goes out of scope and
  // that also triggers deletion of the plugin)
  cs->delete_self_when_unused();

  return result_local;  
}

// check if the jet has explicit_ghosts (knowing that there is an
// area support)
bool Pruner::_check_explicit_ghosts(const PseudoJet &jet) const{
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

// see if there is a common recombiner among the pieces; if there is
// return true and set jet_def_for_recombiner so that the recombiner
// can be taken from that JetDefinition. Otherwise, return
// false. 'assigned' is initially false; when true, each time we meet
// a new jet definition, we'll check it shares the same recombiner as
// jet_def_for_recombiner.
bool Pruner::_check_common_recombiner(const PseudoJet &jet, 
				      JetDefinition &jet_def_for_recombiner,
				      bool assigned) const{
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


// transformer description
std::string Pruner::description() const{
  ostringstream oss;
  oss << "Pruner with jet_definition = (" << _jet_def.description() << ")";
  if (_zcut_dyn) {
    oss << ", dynamic zcut (" << _zcut_dyn->description() << ")"
        << ", dynamic Rcut (" << _Rcut_dyn->description() << ")";
  } else {
    oss << ", zcut = " << _zcut
        << ", Rcut_factor = " << _Rcut_factor;
  }
  return oss.str();
}



//----------------------------------------------------------------------
// class PrunerStructure
//----------------------------------------------------------------------

// return the other jets that may have been found along with the
// result of the pruning
// The resulting vector is sorted in pt
vector<PseudoJet> PrunerStructure::extra_jets() const{ 
  return sorted_by_pt((!SelectorNHardest(1))(validated_cs()->inclusive_jets()));;
}


//----------------------------------------------------------------------
// class PruningRecombiner
//----------------------------------------------------------------------

// decide whether to recombine things or not
void PruningRecombiner::recombine(const PseudoJet &pa, 
                                  const PseudoJet &pb,
                                  PseudoJet &pab) const{
  PseudoJet p;
  _recombiner->recombine(pa, pb, p);

  // if the 2 particles are close enough, do the recombination
  if (pa.squared_distance(pb)<=_Rcut2){
    pab=p; return;
  }

  double pt2a = pa.perp2();
  double pt2b = pb.perp2();

  // check which is the softest
  if (pt2a < pt2b){
    if (pt2a<_zcut2*p.perp2()){
      pab = pb; _rejected.push_back(pa.cluster_hist_index());
    } else {
      pab = p;
    }
  } else {
    if (pt2b<_zcut2*p.perp2()) {
      pab = pa; _rejected.push_back(pb.cluster_hist_index());
    } else {
      pab = p;
    }
  }
}

// description
string PruningRecombiner::description() const{
  ostringstream oss;
  oss << "Pruning recombiner with zcut = " << sqrt(_zcut2)
      << ", Rcut = " << sqrt(_Rcut2)
      << ", and underlying recombiner = " << _recombiner->description();
  return oss.str();
}




//----------------------------------------------------------------------
// class PruningPlugin
//----------------------------------------------------------------------
// the actual clustering work for the plugin
void PruningPlugin::run_clustering(ClusterSequence &input_cs) const{
  // declare a pruning recombiner
  PruningRecombiner pruning_recombiner(_zcut, _Rcut, _jet_def.recombiner());
  JetDefinition jet_def = _jet_def;
  jet_def.set_recombiner(&pruning_recombiner);

  // cluster the particles using that recombiner
  ClusterSequence internal_cs(input_cs.jets(), jet_def);
  const vector<ClusterSequence::history_element> & internal_hist = internal_cs.history();

  // transfer the list of "childless" elements into a bool vector
  vector<bool> kept(internal_hist.size(), true);
  const vector<unsigned int> &pr_rej = pruning_recombiner.rejected();
  for (unsigned int i=0;i<pr_rej.size(); i++) kept[pr_rej[i]]=false;

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

// returns the plugin description
string PruningPlugin::description() const{
  ostringstream oss;
  oss << "Pruning plugin with jet_definition = (" << _jet_def.description()
      <<"), zcut = " << _zcut
      << ", Rcut = " << _Rcut;
  return oss.str();
}


FASTJET_END_NAMESPACE
