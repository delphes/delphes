//FJSTARTHEADER
// $Id: ClusterSequenceActiveAreaExplicitGhosts.cc 3433 2014-07-23 08:17:03Z salam $
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

#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include<limits>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

// save some typing
typedef ClusterSequenceActiveAreaExplicitGhosts ClustSeqActAreaEG;


LimitedWarning ClustSeqActAreaEG::_warnings;

//----------------------------------------------------------------------
///
void ClustSeqActAreaEG::_add_ghosts (
			 const GhostedAreaSpec & ghost_spec) {

  // add the ghosts to the jets
  ghost_spec.add_ghosts(_jets);

  // now add labelling...
  for (unsigned i = _initial_hard_n; i < _jets.size(); i++) {
    //_jets[i].set_user_index(1);
    _is_pure_ghost.push_back(true);
  }

  // and record some info from the ghost_spec
  _ghost_area = ghost_spec.actual_ghost_area();
  _n_ghosts   = ghost_spec.n_ghosts();
}


//----------------------------------------------------------------------
// return the area of a jet
double ClustSeqActAreaEG::area (const PseudoJet & jet) const {
  return _areas[jet.cluster_hist_index()];
}


//----------------------------------------------------------------------
// return the total area
double ClustSeqActAreaEG::total_area () const {
  return _n_ghosts * _ghost_area;
}


//----------------------------------------------------------------------
// return the extended area of a jet
PseudoJet ClustSeqActAreaEG::area_4vector (const PseudoJet & jet) const {
  return _area_4vectors[jet.cluster_hist_index()];
}

//----------------------------------------------------------------------
bool ClustSeqActAreaEG::is_pure_ghost(const PseudoJet & jet) const 
{
  return _is_pure_ghost[jet.cluster_hist_index()];
}

//----------------------------------------------------------------------
bool ClustSeqActAreaEG::is_pure_ghost(int hist_ix) const 
{
  return hist_ix >= 0 ? _is_pure_ghost[hist_ix] : false;
}

//----------------------------------------------------------------------
double ClustSeqActAreaEG::empty_area(const Selector & selector) const {
  // make sure that the selector applies jet by jet
  if (! selector.applies_jet_by_jet()){
    throw Error("ClusterSequenceActiveAreaExplicitGhosts: empty area can only be computed from selectors applying jet by jet");
  }

  vector<PseudoJet> unclust = unclustered_particles();
  double area_local = 0.0;
  for (unsigned iu = 0; iu < unclust.size();  iu++) {
    if (is_pure_ghost(unclust[iu]) && selector.pass(unclust[iu])) {
      area_local += _ghost_area;
    }
  }
  return area_local;
}

//======================================================================
// sort out the areas
void ClustSeqActAreaEG::_post_process() {

  // first check for danger signals.
  // Establish largest ghost transverse momentum
  _max_ghost_perp2 = 0.0;
  for (int i = 0; i < _initial_n; i++) {
    if (_is_pure_ghost[i] && _jets[i].perp2() > _max_ghost_perp2) 
      _max_ghost_perp2 = _jets[i].perp2();
  }

  // now find out if any of the particles are close to danger
  double danger_ratio = numeric_limits<double>::epsilon();
  danger_ratio = danger_ratio * danger_ratio;
  _has_dangerous_particles = false;
  for (int i = 0; i < _initial_n; i++) {
    if (!_is_pure_ghost[i] && 
        danger_ratio * _jets[i].perp2() <=  _max_ghost_perp2) {
      _has_dangerous_particles = true;
      break;
    }
  }

  if (_has_dangerous_particles) _warnings.warn("ClusterSequenceActiveAreaExplicitGhosts: \n  ghosts not sufficiently soft wrt some of the input particles\n  a common cause is (unphysical?) input particles with pt=0 but finite rapidity");

  // sort out sizes
  _areas.resize(_history.size());
  _area_4vectors.resize(_history.size());
  _is_pure_ghost.resize(_history.size());
  
//   copy(_jets.begin(), _jets.begin()+_initial_n, _area_4vectors.begin());
//   for (int i = 0; i < _initial_n; i++) {
//     if (_is_pure_ghost[i]) {
//       _areas[i] = _ghost_area;
//       // normalise pt to be _ghost_area (NB we make use of fact that
//       // for initial particles, jet and clust_hist index are the same).
//       _area_4vectors[i] *= (_ghost_area/_jets[i].perp());
//     } else {
//       _areas[i] = 0;
//       _area_4vectors[i].reset(0,0,0,0);
//     }
//   }
  
  // First set up areas for the initial particles (ghost=_ghost_area,
  // real particles = 0); recall that _initial_n here is the number of
  // particles including ghosts
  for (int i = 0; i < _initial_n; i++) {
    if (_is_pure_ghost[i]) {
      _areas[i] = _ghost_area;
      // normalise pt to be _ghost_area (NB we make use of fact that
      // for initial particles, jet and clust_hist index are the same).
      //_area_4vectors[i] = (_ghost_area/_jets[i].perp()) * _jets[i];
      
      // NB: we use reset_momentum here, to ensure that the area 4
      // vectors do not acquire any structure (the structure would not
      // be meaningful for an area, and it messes up the use count (->
      // memory leaks if the user call delete_self_when_unused).
      _area_4vectors[i].reset_momentum(_jets[i]);
      _area_4vectors[i] *= (_ghost_area/_jets[i].perp());
    } else {
      _areas[i] = 0;
      _area_4vectors[i] = PseudoJet(0.0,0.0,0.0,0.0);
    }
  }
  
  // next follow the branching through and set up the areas 
  // and ghost-nature at each step of the clustering (rather than
  // each jet).
  for (unsigned i = _initial_n; i < _history.size(); i++) {
    if (_history[i].parent2 == BeamJet) {
      _is_pure_ghost[i]  = _is_pure_ghost[_history[i].parent1];
      _areas[i]          = _areas[_history[i].parent1];
      _area_4vectors[i] = _area_4vectors[_history[i].parent1];
    } else {
      _is_pure_ghost[i]  = _is_pure_ghost[_history[i].parent1] && 
	                   _is_pure_ghost[_history[i].parent2]   ;
      _areas[i]          = _areas[_history[i].parent1] + 
	                   _areas[_history[i].parent2]  ;			   
      _jet_def.recombiner()->recombine(_area_4vectors[_history[i].parent1], 
	                   _area_4vectors[_history[i].parent2],
			   _area_4vectors[i]);	
//      _area_4vectors[i] = _area_4vectors[_history[i].parent1] + 
//                          _area_4vectors[_history[i].parent2]  ;
    }

  }
  
}

FASTJET_END_NAMESPACE

