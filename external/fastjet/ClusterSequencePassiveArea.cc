//FJSTARTHEADER
// $Id: ClusterSequencePassiveArea.cc 3433 2014-07-23 08:17:03Z salam $
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

#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceVoronoiArea.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


using namespace std;

//----------------------------------------------------------------------
/// global routine for initialising and running a passive area that is
/// correct in general, but that chooses an optimal approach for
/// various special cases.
void ClusterSequencePassiveArea::_initialise_and_run_PA (
		const JetDefinition & jet_def_in,
		const GhostedAreaSpec & area_spec,
		const bool & writeout_combinations) {

  if (jet_def_in.jet_algorithm() == kt_algorithm) {
    // first run the passive area
    ClusterSequenceVoronoiArea csva(_jets,jet_def_in,VoronoiAreaSpec(1.0));
    // now set up and transfer relevant information    
    // first the clustering sequence
    transfer_from_sequence(csva);
    // then the areas
    _resize_and_zero_AA();
    for (unsigned i = 0; i < _history.size(); i++) {
      int ijetp = _history[i].jetp_index;
      if (ijetp != Invalid) {
        _average_area[i] = csva.area(_jets[ijetp]);
        _average_area_4vector[i] = csva.area_4vector(_jets[ijetp]);
      }
    }

  } else if (jet_def_in.jet_algorithm() == cambridge_algorithm) {
    // run a variant of the cambridge algorithm that has been hacked
    // to deal with passive areas
    JetDefinition tmp_jet_def = jet_def_in;
    tmp_jet_def.set_jet_finder(cambridge_for_passive_algorithm);
    tmp_jet_def.set_extra_param(sqrt(area_spec.mean_ghost_kt()));
    _initialise_and_run_AA(tmp_jet_def, area_spec, writeout_combinations);
    _jet_def = jet_def_in;

  } else if (jet_def_in.jet_algorithm() == antikt_algorithm) {
    // for the antikt algorithm, passive and active are identical
    _initialise_and_run_AA(jet_def_in, area_spec, writeout_combinations);

  } else if (jet_def_in.jet_algorithm() == plugin_algorithm &&
             jet_def_in.plugin()->supports_ghosted_passive_areas()) {
    // for some plugin algorithms, one can "prime" the algorithm with information
    // about the ghost scale, and then an "AA" run will actually give a passive
    // area
    double ghost_sep_scale_store = jet_def_in.plugin()->ghost_separation_scale();
    jet_def_in.plugin()->set_ghost_separation_scale(sqrt(area_spec.mean_ghost_kt()));
    _initialise_and_run_AA(jet_def_in, area_spec, writeout_combinations);

    // restore the original ghost_sep_scale
    jet_def_in.plugin()->set_ghost_separation_scale(ghost_sep_scale_store);

  } else {
    // for a generic algorithm, just run the 1GhostPassiveArea
    _initialise_and_run_1GPA(jet_def_in, area_spec, writeout_combinations);
  }
}

//----------------------------------------------------------------------
// dispatch to most relevant empty area calculation...
double ClusterSequencePassiveArea::empty_area (const Selector & selector) const {
  if (jet_def().jet_algorithm() == kt_algorithm) {
    // run the naive algorithm
    return ClusterSequenceAreaBase::empty_area(selector);
  } else {
    return ClusterSequence1GhostPassiveArea::empty_area(selector);
  }
}


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

