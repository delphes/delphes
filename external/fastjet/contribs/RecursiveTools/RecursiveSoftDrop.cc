// $Id: RecursiveSoftDrop.cc 1192 2018-10-30 16:08:36Z gsoyez $
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

#include "RecursiveSoftDrop.hh"
#include "fastjet/ClusterSequence.hh"

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

namespace internal_recursive_softdrop{
  
//========================================================================
/// \class RSDHistoryElement
/// a helper class to help keeping track od the RSD tree
///
/// The element is created at the top of a branch and updated each
/// time one grooms something away.
class RSDHistoryElement{
public:
  RSDHistoryElement(const PseudoJet &jet, const RecursiveSoftDrop* rsd_ptr, double R0sqr) :
    R0_squared(R0sqr),
    child1_in_history(-1), child2_in_history(-1), symmetry(-1.0), mu2(-1.0){
    reset(jet, rsd_ptr);
  }

  void reset(const PseudoJet &jet, const RecursiveSoftDrop* rsd_ptr){
    current_in_ca_tree = jet.cluster_hist_index();
    PseudoJet piece1, piece2;
    theta_squared = (jet.has_parents(piece1, piece2))
      ? rsd_ptr->squared_geometric_distance(piece1,piece2) : 0.0;
  }
  
  int current_in_ca_tree;  ///< (history) index of the current particle in the C/A tree 
  double theta_squared;    ///< squared angle at which this decays
  double R0_squared;       ///< squared angle at the top of the branch
                           ///< (used for RSD with dynamic_R0)
  int child1_in_history;   ///< hardest of the 2 decay products (-1 if untagged)
  int child2_in_history;   ///< softest of the 2 decay products (-1 if untagged)

  // info about what has been dropped and the local substructure
  vector<double> dropped_delta_R;
  vector<double> dropped_symmetry;
  vector<double> dropped_mu;
  double symmetry, mu2;
};


/// \class OrderRSDHistoryElements
/// angular ordering of (pointers to) the history elements
///
/// our priority queue will use pointers to these elements that are
/// ordered in angle (of the objects they point to)
class OrderRSDHistoryElements{
public:
  bool operator()(const RSDHistoryElement *e1, const RSDHistoryElement *e2) const {
    return e1->theta_squared < e2->theta_squared;
  }
};

} // internal_recursive_softdrop
  
//========================================================================

// initialise all the flags and parameters to their default value
void RecursiveSoftDrop::set_defaults(){
  set_fixed_depth_mode(false);
  set_dynamical_R0(false);
  set_hardest_branch_only(false);
  set_min_deltaR_squared(-1.0);
}

// description of the tool
string RecursiveSoftDrop::description() const{
  ostringstream res;
  res << "recursive application of ["
      << RecursiveSymmetryCutBase::description()
      << "]";
  
  if (_fixed_depth){
    res << ", recursively applied down to a maximal depth of N=";
    if (_n==-1) res << "infinity"; else res << _n;
  } else {
    res << ", applied N=";
    if (_n==-1) res << "infinity"; else res << _n;
    res << " times";
  }
  
  if (_dynamical_R0)
    res << ", with R0 dynamically scaled";
  else
    res << ", with R0 kept fixed";
    
  if (_hardest_branch_only)
    res << ", following only the hardest branch";

  if (_min_dR2>0)
    res << ", with minimal angle (squared) = " << _min_dR2;

  return res.str();
}
  
  
// action on a single jet with RecursiveSoftDrop.
//
// uses "result_fixed_tags" by default (i.e. recurse from R0 to
// smaller angles until n SD conditions have been met), or
// "result_fixed_depth" where each of the previous SD branches are
// recirsed into down to a depth of n.
PseudoJet RecursiveSoftDrop::result(const PseudoJet &jet) const{
  return _fixed_depth ? result_fixed_depth(jet) : result_fixed_tags(jet);
}

// this routine applies the Soft Drop criterion recursively on the
// CA tree until we find n subjets (or until it converges), and
// adds them together into a groomed PseudoJet
PseudoJet RecursiveSoftDrop::result_fixed_tags(const PseudoJet &jet) const {
  // start by reclustering jet with C/A algorithm
  PseudoJet ca_jet = _recluster_if_needed(jet);

  if (! ca_jet.has_valid_cluster_sequence()){
    throw Error("RecursiveSoftDrop can only be applied to jets associated to a (valid) cluster sequence");
  }
  
  const ClusterSequence *cs = ca_jet.validated_cluster_sequence();
  const vector<ClusterSequence::history_element> &cs_history = cs->history();
  const vector<PseudoJet> &cs_jets = cs->jets();

  // initialise counter to 1 subjet (i.e. the full ca_jet)
  int n_tagged = 0;
  int max_njet = ca_jet.constituents().size();

  // create the list of branches
  unsigned int max_history_size = 2*max_njet;
  if ((_n>0) && (_n<max_njet-1)){ max_history_size = 2*(_n+1); }

  // we need to pre-allocate the space for the vector so that the
  // pointers are not invalidated
  vector<internal_recursive_softdrop::RSDHistoryElement> history;
  history.reserve(max_history_size);  // could be one shorter
  history.push_back(internal_recursive_softdrop::RSDHistoryElement(ca_jet, this, _R0sqr));
  
  // create a priority queue containing the subjets and a comparison definition
  priority_queue<internal_recursive_softdrop::RSDHistoryElement*, vector<internal_recursive_softdrop::RSDHistoryElement*>, internal_recursive_softdrop::OrderRSDHistoryElements> active_branches;
  active_branches.push(& (history[0]));

  PseudoJet parent, piece1, piece2;
  double sym, mu2;
  
  // loop over C/A tree until we reach the appropriate number of subjets
  while ((continue_grooming(n_tagged)) && (active_branches.size())) {
    // get the element corresponding to the max dR and the associated PJ
    internal_recursive_softdrop::RSDHistoryElement * elm = active_branches.top();
    PseudoJet parent = cs_jets[cs_history[elm->current_in_ca_tree].jetp_index];
    
    // do one step of SD
    RecursionStatus status = recurse_one_step(parent, piece1, piece2, sym, mu2, &elm->R0_squared);

    // check if we passed the SD condition
    if (status==recursion_success){
      // check for the optional angular cut
      if ((_min_dR2 > 0) && (squared_geometric_distance(piece1,piece2) < _min_dR2))
        break;

      // both subjets are kept in the list for potential further de-clustering
      elm->child1_in_history = history.size();
      elm->child2_in_history = history.size()+1;
      elm->symmetry = sym;
      elm->mu2      = mu2;
      active_branches.pop();

      // update the history
      double next_R0_squared = (_dynamical_R0)
        ? piece1.squared_distance(piece2) : elm->R0_squared;

      internal_recursive_softdrop::RSDHistoryElement elm1(piece1, this, next_R0_squared);
      history.push_back(elm1);
      active_branches.push(&(history.back()));
      internal_recursive_softdrop::RSDHistoryElement elm2(piece2, this, next_R0_squared);
      history.push_back(elm2);
      if (!_hardest_branch_only){
        active_branches.push(&(history.back()));
      }
      
      ++n_tagged;
    } else if (status==recursion_dropped){
       // check for the optional angular cut
      if ((_min_dR2 > 0) && (squared_geometric_distance(piece1,piece2) < _min_dR2))
        break;

      active_branches.pop();
      // tagging failed and the softest branch should be dropped
      // keep track of what has been groomed away
      max_njet -= piece2.constituents().size();
      elm->dropped_delta_R .push_back((elm->theta_squared >= 0) ? sqrt(elm->theta_squared) : -sqrt(elm->theta_squared));
      elm->dropped_symmetry.push_back(sym);
      elm->dropped_mu      .push_back((mu2>=0) ? sqrt(mu2) : -sqrt(mu2));
      
      // keep the hardest branch in the recursion
      elm->reset(piece1, this);
      active_branches.push(elm);
    } else if (status==recursion_no_parents){
      if (_min_dR2 > 0) break;
      active_branches.pop();
      // nothing specific to do: we just keep the curent jet as a "leaf"
    } else { // recursion_issue
      active_branches.pop();
      // we've met an issue
      // if the piece2 is null as well, it means we've had a critical problem.
      // In that case, return an empty PseudoJet
      if (piece2 == 0) return PseudoJet();

      // otherwise, we should consider "piece2" as a final particle
      // not to be recursed into
      if (_min_dR2 > 0) break;
      max_njet -= (piece2.constituents().size()-1);
      break;
    }
    
    // If the missing number of tags is exactly the number of objects
    // we have left in the recursion, stop
    if (n_tagged == max_njet) break;
  }

  // now we have a bunch of history elements that we can use to build the final jet
  vector<PseudoJet> mapped_to_history(history.size());
  unsigned int history_index = history.size();
  do {
    --history_index;
    const internal_recursive_softdrop::RSDHistoryElement & elm = history[history_index];

    // two kinds of events: either just a final leave, potentially with grooming
    // or a brandhing (also with potential grooming at the end)
    if (elm.child1_in_history<0){
      // this is a leaf, i.e. with no further substructure
      PseudoJet & subjet = mapped_to_history[history_index]
        = cs_jets[cs_history[elm.current_in_ca_tree].jetp_index];

      StructureType * structure = new StructureType(subjet);
      if (has_verbose_structure()){
        structure->set_verbose(true);
        structure->set_dropped_delta_R (elm.dropped_delta_R);
        structure->set_dropped_symmetry(elm.dropped_symmetry);
        structure->set_dropped_mu      (elm.dropped_mu);
      }
      subjet.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(structure));
    } else {
      PseudoJet & subjet = mapped_to_history[history_index]
        = join(mapped_to_history[elm.child1_in_history], mapped_to_history[elm.child2_in_history]);
      StructureType * structure = new StructureType(subjet, sqrt(elm.theta_squared), elm.symmetry, sqrt(elm.mu2));
      if (has_verbose_structure()){
        structure->set_verbose(true);
        structure->set_dropped_delta_R (elm.dropped_delta_R);
        structure->set_dropped_symmetry(elm.dropped_symmetry);
        structure->set_dropped_mu      (elm.dropped_mu);
      }
      subjet.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(structure));
    }    
  } while (history_index>0);

  return mapped_to_history[0];
}

// this routine applies the Soft Drop criterion recursively on the
// CA tree, recursing into all the branches found during the previous iteration
// until n layers have been found (or until it converges)
PseudoJet RecursiveSoftDrop::result_fixed_depth(const PseudoJet &jet) const {
  // start by reclustering jet with C/A algorithm
  PseudoJet ca_jet = _recluster_if_needed(jet);

  if (! ca_jet.has_valid_cluster_sequence()){
    throw Error("RecursiveSoftDrop can only be applied to jets associated to a (valid) cluster sequence");
  }
  
  const ClusterSequence *cs = ca_jet.validated_cluster_sequence();
  const vector<ClusterSequence::history_element> &cs_history = cs->history();
  const vector<PseudoJet> &cs_jets = cs->jets();

  // initialise counter to 1 subjet (i.e. the full ca_jet)
  int n_depth = 0;
  int max_njet = ca_jet.constituents().size();

  // create the list of branches
  unsigned int max_history_size = 2*max_njet;
  //if ((_n>0) && (_n<max_njet-1)){ max_history_size = 2*(_n+1); }

  // we need to pre-allocate the space for the vector so that the
  // pointers are not invalidated
  vector<internal_recursive_softdrop::RSDHistoryElement> history;
  history.reserve(max_history_size);  // could be one shorter
  history.push_back(internal_recursive_softdrop::RSDHistoryElement(ca_jet, this, _R0sqr));
  history.back().theta_squared = _R0sqr;
  
  // create a priority queue containing the subjets and a comparison definition
  list<internal_recursive_softdrop::RSDHistoryElement*> active_branches;
  active_branches.push_back(& (history[0]));

  PseudoJet parent, piece1, piece2;
  
  while ((continue_grooming(n_depth)) && (active_branches.size())) {
    // loop over all the branches and look for substructure
    list<internal_recursive_softdrop::RSDHistoryElement*>::iterator hist_it=active_branches.begin();
    while (hist_it!=active_branches.end()){
      // get the element corresponding to the max dR and the associated PJ
      internal_recursive_softdrop::RSDHistoryElement * elm = (*hist_it);
      PseudoJet parent = cs_jets[cs_history[elm->current_in_ca_tree].jetp_index];

      // we need to iterate this branch until we find some substructure
      PseudoJet result_sd;
      if (_dynamical_R0){
        SoftDrop sd(_beta, _symmetry_cut, symmetry_measure(), sqrt(elm->theta_squared),
                    mu_cut(), recursion_choice(), subtractor());
        sd.set_reclustering(false);
        sd.set_verbose_structure(has_verbose_structure());
        result_sd = sd(parent);
      } else {
        result_sd = SoftDrop::result(parent);
      }

      // if we had an empty PJ, that means we ran into some problems.
      // just return an empty PJ ourselves
      if (result_sd == 0) return PseudoJet();
      
      // update the history element to reflect our iteration
      elm->current_in_ca_tree = result_sd.cluster_hist_index();
      
      if (has_verbose_structure()){
        elm->dropped_delta_R    = result_sd.structure_of<SoftDrop>().dropped_delta_R();
        elm->dropped_symmetry   = result_sd.structure_of<SoftDrop>().dropped_symmetry();
        elm->dropped_mu         = result_sd.structure_of<SoftDrop>().dropped_mu();
      }

      // if some substructure was found:
      if (result_sd.structure_of<SoftDrop>().has_substructure()){
        // update the history element to reflect our iteration
        elm->child1_in_history  = history.size();
        elm->child2_in_history  = history.size()+1;
        elm->theta_squared      = result_sd.structure_of<SoftDrop>().delta_R();
        elm->theta_squared     *= elm->theta_squared;
        elm->symmetry           = result_sd.structure_of<SoftDrop>().symmetry();
        elm->mu2                = result_sd.structure_of<SoftDrop>().mu();
        elm->mu2               *= elm->mu2;
             
        // the next iteration will have to handle 2 new history
        // elements (the R0squared argument here is unused)
        result_sd.has_parents(piece1, piece2);
        internal_recursive_softdrop::RSDHistoryElement elm1(piece1, this, _R0sqr);
        history.push_back(elm1);
        // insert it in the active branches if needed
        if (elm1.theta_squared>0)
          active_branches.insert(hist_it,&(history.back())); // insert just before

        internal_recursive_softdrop::RSDHistoryElement elm2(piece2, this, _R0sqr);
        history.push_back(elm2);
        if ((!_hardest_branch_only) && (elm2.theta_squared>0)){
          active_branches.insert(hist_it,&(history.back())); // insert just before
        }
      }
      // otherwise we've just reached the end of the recursion the
      // history information has been updated above
      //
      // we just need to make sure that we do not recurse into that
      // element any longer

      list<internal_recursive_softdrop::RSDHistoryElement*>::iterator current = hist_it;
      ++hist_it;
      active_branches.erase(current);
    } // loop over branches at current depth
    ++n_depth;
  } // loop over depth

  // now we have a bunch of history elements that we can use to build the final jet
  vector<PseudoJet> mapped_to_history(history.size());
  unsigned int history_index = history.size();
  do {
    --history_index;
    const internal_recursive_softdrop::RSDHistoryElement & elm = history[history_index];

    // two kinds of events: either just a final leave, poteitially  with grooming
    // or a brandhing (also with potential grooming at the end)
    if (elm.child1_in_history<0){
      // this is a leaf, i.e. with no further sustructure
      PseudoJet & subjet = mapped_to_history[history_index]
        = cs_jets[cs_history[elm.current_in_ca_tree].jetp_index];

      StructureType * structure = new StructureType(subjet);
      if (has_verbose_structure()){
        structure->set_verbose(true);
      }
      subjet.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(structure));
    } else {
      PseudoJet & subjet = mapped_to_history[history_index]
        = join(mapped_to_history[elm.child1_in_history], mapped_to_history[elm.child2_in_history]);
      StructureType * structure = new StructureType(subjet, sqrt(elm.theta_squared), elm.symmetry, sqrt(elm.mu2));
      if (has_verbose_structure()){
        structure->set_verbose(true);
        structure->set_dropped_delta_R (elm.dropped_delta_R);
        structure->set_dropped_symmetry(elm.dropped_symmetry);
        structure->set_dropped_mu      (elm.dropped_mu);
      }
      subjet.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(structure));
    }    
  } while (history_index>0);

  return mapped_to_history[0];
}


//========================================================================
// implementation of the helpers
//========================================================================

// helper to get all the prongs in a jet that has been obtained using
// RecursiveSoftDrop (instead of recursively parsing the 1->2
// composite jet structure)
vector<PseudoJet> recursive_soft_drop_prongs(const PseudoJet & rsd_jet){
  // make sure that the jet has the appropriate RecursiveSoftDrop structure
  if (!rsd_jet.has_structure_of<RecursiveSoftDrop>())
    return vector<PseudoJet>();

  // if this jet has no substructure, just return a 1-prong object
  if (!rsd_jet.structure_of<RecursiveSoftDrop>().has_substructure())
    return vector<PseudoJet>(1, rsd_jet);
  
  // otherwise fill a vector with all the prongs (no specific ordering)
  vector<PseudoJet> prongs;

  // parse the list of PseudoJet we still need to deal with
  vector<PseudoJet> to_parse = rsd_jet.pieces();  // valid both for a C/A recombination step or a RSD join
  unsigned int i_parse = 0;
  while (i_parse<to_parse.size()){
    const PseudoJet &current = to_parse[i_parse];

    if ((current.has_structure_of<RecursiveSymmetryCutBase>()) &&
        (current.structure_of<RecursiveSymmetryCutBase>().has_substructure())){
      // if this has some deeper substructure, add it to the list of
      // things to further process
      vector<PseudoJet> pieces = current.pieces();
      assert(pieces.size() == 2);
      to_parse[i_parse] = pieces[0];
      to_parse.push_back(pieces[1]);
    } else {
      // no further substructure, just add this as a branch
      prongs.push_back(current);
      ++i_parse;
    }
  }
   
  return prongs;
}
    
}

FASTJET_END_NAMESPACE
