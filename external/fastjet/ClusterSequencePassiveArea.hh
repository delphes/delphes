//FJSTARTHEADER
// $Id: ClusterSequencePassiveArea.hh 3433 2014-07-23 08:17:03Z salam $
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

#ifndef __FASTJET_CLUSTERSEQUENCEPASSIVEAREA_HH__
#define __FASTJET_CLUSTERSEQUENCEPASSIVEAREA_HH__


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence1GhostPassiveArea.hh"
#include<iostream>
#include<vector>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//using namespace std;

/// @ingroup sec_area_classes
/// \class ClusterSequencePassiveArea
/// Like ClusterSequence with computation of the passive jet area
///
/// Class that behaves essentially like ClusterSequence except
/// that it also provides access to the area of a jet (which
/// will be a random quantity... Figure out what to do about seeds 
/// later...)
///
/// This class should not be used directly. Rather use
/// ClusterSequenceArea with the appropriate AreaDefinition
class ClusterSequencePassiveArea : public ClusterSequence1GhostPassiveArea {
public:

  /// constructor based on JetDefinition and PassiveAreaSpec
  template<class L> ClusterSequencePassiveArea
         (const std::vector<L> & pseudojets, 
	  const JetDefinition & jet_def_in,
	  const GhostedAreaSpec & area_spec,
	  const bool & writeout_combinations = false) ;

  /// return an empty area that's appropriate to the passive area
  /// determination carried out
  virtual double empty_area(const Selector & selector) const;

private:

  /// does the initialisation and running specific to the passive
  /// areas class
  void _initialise_and_run_PA (const JetDefinition & jet_def_in,
                               const GhostedAreaSpec & area_spec,
                               const bool & writeout_combinations = false);

};




template<class L> ClusterSequencePassiveArea::ClusterSequencePassiveArea 
(const std::vector<L> & pseudojets, 
 const JetDefinition & jet_def_in,
 const GhostedAreaSpec & area_spec,
 const bool & writeout_combinations) {

  // transfer the initial jets (type L) into our own array
  _transfer_input_jets(pseudojets);

  // run the clustering for passive areas
  _initialise_and_run_PA(jet_def_in, area_spec, writeout_combinations);

}


  
FASTJET_END_NAMESPACE

#endif // __FASTJET_CLUSTERSEQUENCEPASSIVEAREA_HH__
