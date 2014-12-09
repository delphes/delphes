// $Id: Recluster.cc 3629 2014-08-14 17:21:15Z salam $
//
// Copyright (c) 2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#include "fastjet/tools/Recluster.hh"
#include "fastjet/CompositeJetStructure.hh"
#include <fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh>
#include <sstream>
#include <typeinfo>

using namespace std;

// Comments:
//
//  - If the jet comes from a C/A clustering (or is a composite jet
//    made of C/A clusterings) and we recluster it with a C/A
//    algorithm, we just use exclusive jets instead of doing the
//    clustering explicitly. In this specific case, we need to make
//    sure that all the pieces share the same cluster sequence.
//
//  - If the recombiner has to be taken from the original jet and that
//    jet is composite, we need to check that all the pieces were
//    obtained with the same recombiner.
//
//  - Note that a preliminary version of this code has been
//    distributed in the RecursiveTools fastjet-contrib. The present
//    version has a few minor fixed

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

LimitedWarning Recluster::_explicit_ghost_warning;

// ctor
Recluster::Recluster(JetAlgorithm new_jet_alg, double new_jet_radius, Keep keep_in)
  : _new_jet_def(JetDefinition(new_jet_alg, new_jet_radius)), 
    _acquire_recombiner(true), _keep(keep_in), _cambridge_optimisation_enabled(true){}

Recluster::Recluster(JetAlgorithm new_jet_alg, Keep keep_in)
    : _acquire_recombiner(true), _keep(keep_in), _cambridge_optimisation_enabled(true){
  switch (JetDefinition::n_parameters_for_algorithm(new_jet_alg)){
  case 0: _new_jet_def = JetDefinition(new_jet_alg); break;
  case 1: _new_jet_def = JetDefinition(new_jet_alg, JetDefinition::max_allowable_R); break;
  default:
    throw Error("Recluster(): tried to construct specifying only a jet algorithm ("+JetDefinition::algorithm_description(new_jet_alg)+") which takes more than 1 parameter");
  };
}


// class description
string Recluster::description() const {
  ostringstream ostr;
  ostr << "Recluster with new_jet_def = ";
  if (_acquire_recombiner){
    ostr << _new_jet_def.description_no_recombiner();
    ostr << ", using a recombiner obtained from the jet being reclustered";
  } else {
    ostr << _new_jet_def.description();
  }

  if (_keep == keep_only_hardest)
    ostr << " and keeping the hardest inclusive jet";
  else
    ostr << " and joining all inclusive jets into a composite jet";

  return ostr.str();
}


// the main piece of code that performs the reclustering
PseudoJet Recluster::result(const PseudoJet &jet) const {
  // get the incljets and the exact jet definition that has been used
  // to get them
  vector<PseudoJet> incljets;
  bool ca_optimised = get_new_jets_and_def(jet, incljets);

  return generate_output_jet(incljets, ca_optimised);
}


// a lower-level method that does the actual work of reclustering the
// input jet. The resulting incljets are stored in output_jets and the
// jet definition that has been used can be deduced from their
// associated ClusterSequence
//
//  - input_jet       the (input) jet that one wants to recluster
//  - output_jets     incljets resulting from the new clustering
//
// returns true if the C/A optimisation has been used (this means
// that generate_output_jet will watch out for non-explicit-ghost
// areas that might be leftover)
bool Recluster::get_new_jets_and_def(const PseudoJet & input_jet, 
                                     vector<PseudoJet> & output_jets) const{
  // generic sanity checks
  //-------------------------------------------------------------------
  // make sure that the jet has constituents
  if (! input_jet.has_constituents())
    throw Error("Recluster can only be applied on jets having constituents");

  // tests particular to certain configurations
  //-------------------------------------------------------------------

  // for a whole variety of tests, we shall need the "recursive"
  // pieces of the jet (the jet itself or recursing down to its most
  // fundamental pieces). So we do compute these once and for all.
  //
  // Note that the pieces are always needed (either for C/A or for the
  // area checks)
  vector<PseudoJet> all_pieces;
  if ((!_get_all_pieces(input_jet, all_pieces)) || (all_pieces.size()==0)){
    throw Error("Recluster: failed to retrieve all the pieces composing the jet.");
  }

  // decide which jet definition to use
  //-------------------------------------------------------------------
  JetDefinition new_jet_def = _new_jet_def;
  if (_acquire_recombiner){
    _acquire_recombiner_from_pieces(all_pieces, new_jet_def);
  }

  // the vector that will ultimately hold the incljets
  output_jets.clear();

  // check if we can apply the simplification for C/A jets reclustered
  // with C/A
  //
  // we apply C/A clustering iff
  //  - the requested new_jet_def is C/A
  //  - the jet is either directly coming from C/A or if it is a
  //    superposition of C/A jets from the same cluster sequence
  //  - the pieces agree with the recombination scheme of new_jet_def
  //
  // In this case area support will be automatically inherited so we
  // can only worry about this later
  // -------------------------------------------------------------------
  if (_check_ca(all_pieces, new_jet_def)){
    _recluster_ca(all_pieces, output_jets, new_jet_def.R());
    output_jets = sorted_by_pt(output_jets);
    return true;
  }

  // decide if area support has to be kept
  //-------------------------------------------------------------------
  bool include_area_support = input_jet.has_area();
  if ((include_area_support) &&  (!_check_explicit_ghosts(all_pieces))){
    _explicit_ghost_warning.warn("Recluster: the original cluster sequence is lacking explicit ghosts; area support will no longer be available after re-clustering");
    include_area_support = false;
  }

  // extract the incljets
  //-------------------------------------------------------------------
  _recluster_generic(input_jet, output_jets, new_jet_def, include_area_support);
  output_jets = sorted_by_pt(output_jets);

  return false;
}

// given a set of incljets and a jet definition used, create the
// resulting PseudoJet
PseudoJet Recluster::generate_output_jet(std::vector<PseudoJet> & incljets, 
                                         bool ca_optimisation_used) const{
  // first handle the case where we only need to keep the hardest incljet
  if (_keep == keep_only_hardest) {
    if (incljets.size() > 0) {
      return incljets[0];
    } else {
      return PseudoJet();
    }
  }

  // now the case where all incljets have to be joined

  // safekeeper
  if (incljets.size()==0) return join(incljets);

  PseudoJet reclustered = join(incljets, 
			       *(incljets[0].associated_cluster_sequence()->jet_def().recombiner()));

  // if we've used C/A optimisation, we need to get rid of the area
  // information if it comes from a non-explicit-ghost clustering.
  // (because in that case it can be erroneous due to the lack of
  // information about empty areas)
  if (ca_optimisation_used){
    if (reclustered.has_area() &&
        (incljets.size() > 0) &&
        (! incljets[0].validated_csab()->has_explicit_ghosts())){
      CompositeJetStructure *css = dynamic_cast<CompositeJetStructure *>(reclustered.structure_non_const_ptr());
      assert(css);
      css->discard_area();
    }
  }
  
  return reclustered;
}

//----------------------------------------------------------------------
// the parts that really do the reclustering
//----------------------------------------------------------------------

// get the subjets in the simple case of C/A+C/A
void Recluster::_recluster_ca(const vector<PseudoJet> & all_pieces, 
                              vector<PseudoJet> & subjets, 
                              double Rfilt) const{
  subjets.clear();

  // each individual piece should have a C/A cluster sequence
  // associated to it. So a simple loop would do the job
  for (vector<PseudoJet>::const_iterator piece_it = all_pieces.begin(); 
       piece_it!=all_pieces.end(); piece_it++){
    // just extract the exclusive subjets of 'jet'
    const ClusterSequence *cs = piece_it->associated_cluster_sequence(); 
    vector<PseudoJet> local_subjets;

    double dcut = Rfilt / cs->jet_def().R();
    if (dcut>=1.0){
      // remember that in this case all the pairwise interpiece
      // distances are supposed to be larger than Rfilt (this was
      // tested in _check_ca), which means that they can never
      // recombine with each other.
      local_subjets.push_back(*piece_it);
    } else {
      local_subjets = piece_it->exclusive_subjets(dcut*dcut);
    }

    copy(local_subjets.begin(), local_subjets.end(), back_inserter(subjets));
  }
}


// perform the reclustering itself for all cases where the "C/A trick"
// does not apply
void Recluster::_recluster_generic(const PseudoJet & jet, 
                                   vector<PseudoJet> & incljets,
                                   const JetDefinition & new_jet_def,
                                   bool do_areas) const{
  // create a new, internal, ClusterSequence from the jet constituents
  // get the incljets directly from there
  //
  // If the jet has area support then we separate the ghosts from the
  // "regular" particles so the incljets will also have area
  // support. Note that we do this regardless of whether rho is zero
  // or not.
  //
  // Note that to be able to separate the ghosts, one needs explicit
  // ghosts!!
  // ---------------------------------------------------------------
  if (do_areas){
    //vector<PseudoJet> all_constituents = jet.constituents();
    vector<PseudoJet> regular_constituents, ghosts;  
    SelectorIsPureGhost().sift(jet.constituents(), ghosts, regular_constituents);

    // figure out the ghost area from the 1st ghost (if none, any value
    // would probably do as the area will be 0 and subtraction will have
    // no effect!)
    double ghost_area = (ghosts.size()) ? ghosts[0].area() : 0.01;
    ClusterSequenceActiveAreaExplicitGhosts * csa
      = new ClusterSequenceActiveAreaExplicitGhosts(regular_constituents, 
                                                    new_jet_def, 
                                                    ghosts, ghost_area);

    incljets = csa->inclusive_jets();

    // allow the cs to be deleted when it's no longer used
    // 
    // Note that there is at least one constituent in the jet so there
    // is in principle at least one incljet. But one may have used a
    // nasty recombiner or jet def that left an empty set of incljets,
    // so we'd rather play it safe (e.g. GridJetPlugin, with the
    // constituents outside the range of the grid)
    if (incljets.size())
      csa->delete_self_when_unused();
    else
      delete csa;
  } else {
    ClusterSequence * cs = new ClusterSequence(jet.constituents(), new_jet_def);
    incljets = cs->inclusive_jets();
    // allow the cs to be deleted when it's no longer used (again, we
    // add an extra safety check)
    if (incljets.size())
      cs->delete_self_when_unused();
    else 
      delete cs;
  }
}

//----------------------------------------------------------------------
// various checks and internal constructs
//----------------------------------------------------------------------

// fundamental info for CompositeJets
//----------------------------------------------------------------------
 
// get the pieces down to the fundamental pieces
// 
// Note that this just checks that there is an associated CS to the
// fundamental pieces, not that it is still valid
bool Recluster::_get_all_pieces(const PseudoJet &jet, vector<PseudoJet> &all_pieces) const{
  if (jet.has_associated_cluster_sequence()){
    all_pieces.push_back(jet);
    return true;
  }

  if (jet.has_pieces()){
    const vector<PseudoJet> pieces = jet.pieces();
    for (vector<PseudoJet>::const_iterator it=pieces.begin(); it!=pieces.end(); it++)
      if (!_get_all_pieces(*it, all_pieces)) return false;
    return true;
  }

  return false;
}

// construct the re-clustering jet definition using the recombiner
// from whatever definition has been used to obtain the original jet
//----------------------------------------------------------------------
void Recluster::_acquire_recombiner_from_pieces(const vector<PseudoJet> &all_pieces, 
                                                JetDefinition &new_jet_def) const{
  // a quick safety check
  assert(_acquire_recombiner);

  // check that all the pieces have the same recombiner
  //
  // Note that if the jet has an associated cluster sequence that is no
  // longer valid, an error will be thrown (needed since it could be the
  // 1st check called after the enumeration of the pieces)
  const JetDefinition & jd_ref = all_pieces[0].validated_cs()->jet_def();
  for (unsigned int i=1; i<all_pieces.size(); i++){
    if (!all_pieces[i].validated_cs()->jet_def().has_same_recombiner(jd_ref)){
      throw Error("Recluster instance is configured to determine the recombination scheme (or recombiner) from the original jet, but different pieces of the jet were found to have non-equivalent recombiners.");
    }
  }

  // get the recombiner from the original jet_def
  new_jet_def.set_recombiner(jd_ref);
}

// area support
//----------------------------------------------------------------------

// check if the jet (or all its pieces) have explicit ghosts
// (assuming the jet has area support).
//
// Note that if the jet has an associated cluster sequence that is no
// longer valid, an error will be thrown (needed since it could be the
// 1st check called after the enumeration of the pieces)
bool Recluster::_check_explicit_ghosts(const vector<PseudoJet> &all_pieces) const{
  for (vector<PseudoJet>::const_iterator it=all_pieces.begin(); it!=all_pieces.end(); it++)
    if (! it->validated_csab()->has_explicit_ghosts()) return false;
  return true;
}

// C/A specific tests
//----------------------------------------------------------------------

// check if one can apply the simplification for C/A incljets
//
// This includes:
//  - the incljet definition asks for C/A incljets
//  - all the pieces share the same CS
//  - that CS is C/A with the same recombiner as the incljet def
//  - the re-clustering radius is not larger than any of the pairwise
//    distance between the pieces
//
// Note that if the jet has an associated cluster sequence that is no
// longer valid, an error will be thrown (needed since it could be the
// 1st check called after the enumeration of the pieces)
bool Recluster::_check_ca(const vector<PseudoJet> &all_pieces, 
                          const JetDefinition &new_jet_def) const{
  // check that optimisation is enabled
  if (!_cambridge_optimisation_enabled) return false;

  // check that we're reclustering with C/A
  if (new_jet_def.jet_algorithm() != cambridge_algorithm) return false;

  // check that the 1st of all the pieces (we're sure there is at
  // least one) is coming from a C/A clustering. Then check that all
  // the following pieces share the same ClusterSequence
  const ClusterSequence * cs_ref = all_pieces[0].validated_cs();
  if (cs_ref->jet_def().jet_algorithm() != cambridge_algorithm) return false;
  for (unsigned int i=1; i<all_pieces.size(); i++)
    if (all_pieces[i].validated_cs() != cs_ref) return false;

  // check that the 1st piece has the same recombiner as the one used
  // for the incljet clustering
  // Note that since they share the same CS, checking the 1st one is enough
  if (!cs_ref->jet_def().has_same_recombiner(new_jet_def)) return false;

  // we also have to make sure that the reclustering radius is not larger
  // than any of the inter-piece distances
  double Rnew2 = new_jet_def.R();
  Rnew2 *= Rnew2;
  for (unsigned int i=0; i<all_pieces.size()-1; i++){
    for (unsigned int j=i+1; j<all_pieces.size(); j++){
      if (all_pieces[i].squared_distance(all_pieces[j]) <  Rnew2) return false;
    }
  }

  return true;
}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
