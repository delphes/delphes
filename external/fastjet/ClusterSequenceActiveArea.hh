//FJSTARTHEADER
// $Id: ClusterSequenceActiveArea.hh 3619 2014-08-13 14:17:19Z salam $
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

#ifndef __FASTJET_CLUSTERSEQUENCEACTIVEAREA_HH__
#define __FASTJET_CLUSTERSEQUENCEACTIVEAREA_HH__


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include<iostream>
#include<vector>

//------------ backwards compatibility with version 2.1 -------------
// for backwards compatibility make ActiveAreaSpec name available
//#include "fastjet/ActiveAreaSpec.hh"
//#include "fastjet/ClusterSequenceWithArea.hh"
//--------------------------------------------------------------------


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//using namespace std;

/// @ingroup sec_area_classes
/// \class ClusterSequenceActiveArea
/// Like ClusterSequence with computation of the active jet area
///
/// Class that behaves essentially like ClusterSequence except
/// that it also provides access to the area of a jet (which
/// will be a random quantity... Figure out what to do about seeds 
/// later...)
///
/// This class should not be used directly. Rather use
/// ClusterSequenceArea with the appropriate AreaDefinition
class ClusterSequenceActiveArea : public ClusterSequenceAreaBase {
public:

  /// default constructor
  ClusterSequenceActiveArea() {}

  /// constructor based on JetDefinition and GhostedAreaSpec
  template<class L> ClusterSequenceActiveArea
         (const std::vector<L> & pseudojets, 
	  const JetDefinition & jet_def_in,
	  const GhostedAreaSpec & ghost_spec,
	  const bool & writeout_combinations = false) ;

  virtual double area (const PseudoJet & jet) const {
                             return _average_area[jet.cluster_hist_index()];};
  virtual double area_error (const PseudoJet & jet) const {
                             return _average_area2[jet.cluster_hist_index()];};

  virtual PseudoJet area_4vector (const PseudoJet & jet) const {
                    return _average_area_4vector[jet.cluster_hist_index()];};

  /// enum providing a variety of tentative strategies for estimating
  /// the background (e.g. non-jet) activity in a highly populated event; the
  /// one that has been most extensively tested is median.
  /// 
  /// These strategies are OBSOLETE and deprecated (see comment
  /// for pt_per_unit_area).
  enum mean_pt_strategies{median=0, non_ghost_median, pttot_over_areatot, 
			  pttot_over_areatot_cut, mean_ratio_cut, play,
			  median_4vector};

  /// return the transverse momentum per unit area according to one
  /// of the above strategies; for some strategies (those with "cut"
  /// in their name) the parameter "range" allows one to exclude a
  /// subset of the jets for the background estimation, those that
  /// have pt/area > median(pt/area)*range.
  ///
  /// NB: This call is OBSOLETE and deprecated; use a
  /// JetMedianBackgroundEstimator or GridMedianBackgroundEstimator
  /// instead.
  double pt_per_unit_area(mean_pt_strategies strat=median, 
                          double range=2.0 ) const;

  /// rewrite the empty area from the parent class, so as to use
  /// all info at our disposal
  /// return the total area, corresponding to a given Selector, that
  /// consists of ghost jets or unclustered ghosts
  ///
  /// The selector passed as an argument needs to apply jet by jet.
  virtual double empty_area(const Selector & selector) const;

  /// return the true number of empty jets (replaces
  /// ClusterSequenceAreaBase::n_empty_jets(...))
  virtual double n_empty_jets(const Selector & selector) const;

protected:
  void _resize_and_zero_AA ();
  void _initialise_AA(const JetDefinition & jet_def,
                      const GhostedAreaSpec & ghost_spec,
                      const bool & writeout_combinations,
                      bool & continue_running);

  void _run_AA(const GhostedAreaSpec & ghost_spec);

  void _postprocess_AA(const GhostedAreaSpec & ghost_spec);

  /// does the initialisation and running specific to the active
  /// areas class
  void _initialise_and_run_AA (const JetDefinition & jet_def,
                               const GhostedAreaSpec & ghost_spec,
                               const bool & writeout_combinations = false);

  /// transfer the history (and jet-momenta) from clust_seq to our
  /// own internal structure while removing ghosts
  void _transfer_ghost_free_history(
           const ClusterSequenceActiveAreaExplicitGhosts & clust_seq);


  /// transfer areas from the ClusterSequenceActiveAreaExplicitGhosts
  /// object into our internal area bookkeeping...
  void _transfer_areas(const std::vector<int> & unique_hist_order, 
                       const ClusterSequenceActiveAreaExplicitGhosts & );

  /// child classes benefit from having these at their disposal
  std::valarray<double> _average_area, _average_area2;
  std::valarray<PseudoJet> _average_area_4vector;

  /// returns true if there are any particles whose transverse momentum
  /// if so low that there's a risk of the ghosts having modified the
  /// clustering sequence
  bool has_dangerous_particles() const {return _has_dangerous_particles;}

private:


  double           _non_jet_area, _non_jet_area2, _non_jet_number;

  double _maxrap_for_area; // max rap where we put ghosts
  double _safe_rap_for_area; // max rap where we trust jet areas

  bool   _has_dangerous_particles; 


  /// routine for extracting the tree in an order that will be independent
  /// of any degeneracies in the recombination sequence that don't
  /// affect the composition of the final jets
  void _extract_tree(std::vector<int> &) const;
  /// do the part of the extraction associated with pos, working
  /// through its children and their parents
  void _extract_tree_children(int pos, std::valarray<bool> &, const std::valarray<int> &, std::vector<int> &) const;
  /// do the part of the extraction associated with the parents of pos.
  void _extract_tree_parents (int pos, std::valarray<bool> &, const std::valarray<int> &,  std::vector<int> &) const;

  /// check if two jets have the same momentum to within the
  /// tolerance (and if pt's are not the same we're forgiving and
  /// look to see if the energy is the same)
  void _throw_unless_jets_have_same_perp_or_E(const PseudoJet & jet, 
                                              const PseudoJet & refjet, 
                                              double tolerance,
             const ClusterSequenceActiveAreaExplicitGhosts & jets_ghosted_seq
                                              ) const;

  /// since we are playing nasty games with seeds, we should warn
  /// the user a few times
  //static int _n_seed_warnings;
  //const static int _max_seed_warnings = 10;

  // record the number of repeats
  int _ghost_spec_repeat;

  /// a class for our internal storage of ghost jets
  class GhostJet : public PseudoJet {
  public:
    GhostJet(const PseudoJet & j, double a) : PseudoJet(j), area(a){}
    double area;
  };

  std::vector<GhostJet> _ghost_jets;
  std::vector<GhostJet> _unclustered_ghosts;
};




template<class L> ClusterSequenceActiveArea::ClusterSequenceActiveArea 
(const std::vector<L> & pseudojets, 
 const JetDefinition & jet_def_in,
 const GhostedAreaSpec & ghost_spec,
 const bool & writeout_combinations) {

  // transfer the initial jets (type L) into our own array
  _transfer_input_jets(pseudojets);

  // run the clustering for active areas
  _initialise_and_run_AA(jet_def_in, ghost_spec, writeout_combinations);

}


  
FASTJET_END_NAMESPACE

#endif // __FASTJET_CLUSTERSEQUENCEACTIVEAREA_HH__
