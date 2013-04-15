//STARTHEADER
// $Id: ClusterSequenceActiveAreaExplicitGhosts.hh 2687 2011-11-14 11:17:51Z soyez $
//
// Copyright (c) 2005-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER

#ifndef __FASTJET_CLUSTERSEQUENCEACTIVEAREAEXPLICITGHOSTS_HH_
#define __FASTJET_CLUSTERSEQUENCEACTIVEAREAEXPLICITGHOSTS_HH_ 

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/LimitedWarning.hh"
#include<iostream>
#include<vector>
#include <cstdio>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//======================================================================
/// @ingroup sec_area_classes
/// \class ClusterSequenceActiveAreaExplicitGhosts
/// Like ClusterSequence with computation of the active jet area with the
/// addition of explicit ghosts
///
/// Class that behaves essentially like ClusterSequence except
/// that it also provides access to the area of a jet (which
/// will be a random quantity... Figure out what to do about seeds 
/// later...)
///
/// This class should not be used directly. Rather use
/// ClusterSequenceArea with the appropriate AreaDefinition
class ClusterSequenceActiveAreaExplicitGhosts : 
  public ClusterSequenceAreaBase {
public:
  /// constructor using a GhostedAreaSpec to specify how the area is
  /// to be measured
  template<class L> ClusterSequenceActiveAreaExplicitGhosts
         (const std::vector<L> & pseudojets, 
          const JetDefinition & jet_def_in,
	  const GhostedAreaSpec & ghost_spec,
	  const bool & writeout_combinations = false) 
	   : ClusterSequenceAreaBase() {
           std::vector<L> * ghosts = NULL;
	   _initialise(pseudojets,jet_def_in,&ghost_spec,ghosts,0.0,
                       writeout_combinations); }

  template<class L> ClusterSequenceActiveAreaExplicitGhosts
         (const std::vector<L> & pseudojets, 
          const JetDefinition & jet_def_in,
          const std::vector<L> & ghosts,
          double ghost_area,
	  const bool & writeout_combinations = false) 
	   : ClusterSequenceAreaBase() {
           const GhostedAreaSpec * ghost_spec = NULL;
	   _initialise(pseudojets,jet_def_in,ghost_spec,&ghosts,ghost_area,
                       writeout_combinations); }


  /// does the actual work of initialisation
  template<class L> void _initialise
         (const std::vector<L> & pseudojets, 
          const JetDefinition & jet_def_in,
	  const GhostedAreaSpec * ghost_spec,
	  const std::vector<L> * ghosts,
	  double                 ghost_area,
	  const bool & writeout_combinations); 

  //vector<PseudoJet> constituents (const PseudoJet & jet) const;

  /// returns the number of hard particles (i.e. those supplied by the user).
  unsigned int n_hard_particles() const;

  /// returns the area of a jet
  virtual double area (const PseudoJet & jet) const;

  /// returns a four vector corresponding to the sum (E-scheme) of the
  /// ghost four-vectors composing the jet area, normalised such that
  /// for a small contiguous area the p_t of the extended_area jet is
  /// equal to area of the jet.
  virtual PseudoJet area_4vector (const PseudoJet & jet) const;

  /// true if a jet is made exclusively of ghosts
  virtual bool is_pure_ghost(const PseudoJet & jet) const;

  /// true if the entry in the history index corresponds to a
  /// ghost; if hist_ix does not correspond to an actual particle
  /// (i.e. hist_ix < 0), then the result is false.
  bool is_pure_ghost(int history_index) const;

  /// this class does have explicit ghosts
  virtual bool has_explicit_ghosts() const {return true;}

  /// return the total area, corresponding to a given Selector, that
  /// consists of unclustered ghosts
  ///
  /// The selector needs to apply jet by jet
  virtual double empty_area(const Selector & selector) const;

  /// returns the total area under study
  double total_area () const;
  
  /// returns the largest squared transverse momentum among
  /// all ghosts
  double max_ghost_perp2() const {return _max_ghost_perp2;}

  /// returns true if there are any particles whose transverse momentum
  /// if so low that there's a risk of the ghosts having modified the
  /// clustering sequence
  bool has_dangerous_particles() const {return _has_dangerous_particles;}

  /// get the area of the ghosts
  //double ghost_area() const{return _ghost_area;}

private:

  int    _n_ghosts;
  double _ghost_area;
  std::vector<bool> _is_pure_ghost;
  std::vector<double> _areas;
  std::vector<PseudoJet> _area_4vectors;
  
  // things related to checks for dangerous particles
  double _max_ghost_perp2;
  bool   _has_dangerous_particles; 
  static LimitedWarning _warnings;

  //static int _n_warn_dangerous_particles;
  //static const int _max_warn_dangerous_particles = 5;

  
  unsigned int _initial_hard_n;

  /// adds the "ghost" momenta, which will be used to estimate
  /// the jet area
  void _add_ghosts(const GhostedAreaSpec & ghost_spec); 

  /// another way of adding ghosts
  template<class L> void _add_ghosts (
	  const std::vector<L> & ghosts,
	  double                 ghost_area);

  /// routine to be called after the processing is done so as to
  /// establish summary information on all the jets (areas, whether
  /// pure ghost, etc.)
  void _post_process();

};


//----------------------------------------------------------------------
// initialise from some generic type... Has to be made available
// here in order for the template aspect of it to work...
template<class L> void ClusterSequenceActiveAreaExplicitGhosts::_initialise
         (const std::vector<L> & pseudojets, 
          const JetDefinition & jet_def_in,
	  const GhostedAreaSpec * ghost_spec,
	  const std::vector<L> * ghosts,
	  double                 ghost_area,
	  const bool & writeout_combinations) {
  // don't reserve space yet -- will be done below

  // insert initial jets this way so that any type L that can be
  // converted to a pseudojet will work fine (basically PseudoJet
  // and any type that has [] subscript access to the momentum
  // components, such as CLHEP HepLorentzVector).
  for (unsigned int i = 0; i < pseudojets.size(); i++) {
    PseudoJet mom(pseudojets[i]);
    //mom.set_user_index(0); // for user's particles (user index now lost...)
    _jets.push_back(mom);
    _is_pure_ghost.push_back(false);
  }

  _initial_hard_n = _jets.size();

  if (ghost_spec != NULL) {
    //std::cout << "about to reserve " << (_jets.size()+ghost_spec->n_ghosts())*2 << std::endl;
    _jets.reserve((_jets.size()+ghost_spec->n_ghosts()));
    _add_ghosts(*ghost_spec);
  } else {
    _jets.reserve(_jets.size()+ghosts->size());
    _add_ghosts(*ghosts, ghost_area);
  }

  if (writeout_combinations) {
    std::cout << "# Printing particles including ghosts\n";
    for (unsigned j = 0; j < _jets.size(); j++) {
      printf("%5u %20.13f %20.13f %20.13e\n",
	       j,_jets[j].rap(),_jets[j].phi_02pi(),_jets[j].kt2());
    }
    std::cout << "# Finished printing particles including ghosts\n";
  }

  // this will ensure that we can still point to jets without
  // difficulties arising!
  //std::cout << _jets.size() << " " << _jets.size()*2 << " " << _jets.max_size() << std::endl;
  _jets.reserve(_jets.size()*2); //GPS tmp removed

  // run the clustering
  _initialise_and_run(jet_def_in,writeout_combinations);

  // set up all other information
  _post_process();
}


inline unsigned int ClusterSequenceActiveAreaExplicitGhosts::n_hard_particles() const {return _initial_hard_n;}


//----------------------------------------------------------------------
/// add an explicitly specified bunch of ghosts
template<class L> void ClusterSequenceActiveAreaExplicitGhosts::_add_ghosts (
	  const std::vector<L> & ghosts,
	  double                 ghost_area) {

  
  for (unsigned i = 0; i < ghosts.size(); i++) {
    _is_pure_ghost.push_back(true);
    _jets.push_back(ghosts[i]);
  }
  // and record some info about ghosts
  _ghost_area = ghost_area;
  _n_ghosts   = ghosts.size();
}


FASTJET_END_NAMESPACE

#endif // __FASTJET_CLUSTERSEQUENCEACTIVEAREAEXPLICITGHOSTS_HH_
