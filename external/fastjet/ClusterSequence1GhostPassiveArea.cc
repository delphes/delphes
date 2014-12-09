//FJSTARTHEADER
// $Id: ClusterSequence1GhostPassiveArea.cc 3433 2014-07-23 08:17:03Z salam $
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

#include "fastjet/ClusterSequence1GhostPassiveArea.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


using namespace std;

//----------------------------------------------------------------------
/// global routine for initialising and running a general passive area
void ClusterSequence1GhostPassiveArea::_initialise_and_run_1GPA (
		const JetDefinition & jet_def_in,
		const GhostedAreaSpec & area_spec,
		const bool & writeout_combinations) {

  bool continue_running;
  _initialise_AA(jet_def_in, area_spec, writeout_combinations, continue_running);
  if (continue_running) {
    _run_1GPA(area_spec);
    _postprocess_AA(area_spec);
  }
}


//----------------------------------------------------------------------
/// routine for running a passive area one ghost at a time.
void ClusterSequence1GhostPassiveArea::_run_1GPA (const GhostedAreaSpec & area_spec) {
    // record the input jets as they are currently
  vector<PseudoJet> input_jets(_jets);

  // code for testing the unique tree
  vector<int> unique_tree;

  // initialise our temporary average area^2
  valarray<double> lcl_average_area2(0.0, _average_area.size());
  valarray<double> last_average_area(0.0, _average_area.size());

  // run the clustering multiple times so as to get areas of all the jets
  for (int irepeat = 0; irepeat < area_spec.repeat(); irepeat++) {

    // first establish list of all ghosts
    vector<PseudoJet> all_ghosts;
    area_spec.add_ghosts(all_ghosts);

    // then run many subsets of ghosts (actually each subset contains just one ghost)
    for (unsigned ig = 0; ig < all_ghosts.size(); ig++) {
      vector<PseudoJet> some_ghosts;
      some_ghosts.push_back(all_ghosts[ig]);
      ClusterSequenceActiveAreaExplicitGhosts clust_seq(input_jets, jet_def(), 
                                                       some_ghosts, area_spec.actual_ghost_area());

      if (irepeat == 0 && ig == 0) {
        // take the non-ghost part of the history and put into our own
        // history.
        _transfer_ghost_free_history(clust_seq);
        // get the "unique" order that will be used for transferring all areas. 
        unique_tree = unique_history_order();
      }
      
      // transfer areas from clust_seq into our object
      _transfer_areas(unique_tree, clust_seq);
    }
    // keep track of fluctuations in area
    lcl_average_area2 += (_average_area - last_average_area)*
                         (_average_area - last_average_area);
    last_average_area = _average_area;
  }
  _average_area2 = lcl_average_area2;
}


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

