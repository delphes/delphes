//STARTHEADER
// $Id: ClusterSequence1GhostPassiveArea.hh 2687 2011-11-14 11:17:51Z soyez $
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

#ifndef __FASTJET_CLUSTERSEQUENCE1GHOSTPASSIVEAREA_HH__
#define __FASTJET_CLUSTERSEQUENCE1GHOSTPASSIVEAREA_HH__


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include<iostream>
#include<vector>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//using namespace std;

/// @ingroup sec_area_classes
/// \class ClusterSequence1GhostPassiveArea
/// Like ClusterSequence with computation of the passive jet area by
/// adding a single ghost
///
/// Class that behaves essentially like ClusterSequence except
/// that it also provides access to the area of a jet (which
/// will be a random quantity... Figure out what to do about seeds 
/// later...)
///
/// This class should not be used directly. Rather use
/// ClusterSequenceArea
class ClusterSequence1GhostPassiveArea : public ClusterSequenceActiveArea {
public:

  ClusterSequence1GhostPassiveArea() {}

  /// constructor based on JetDefinition and 1GhostPassiveAreaSpec
  template<class L> ClusterSequence1GhostPassiveArea
         (const std::vector<L> & pseudojets, 
	  const JetDefinition & jet_def_in,
	  const GhostedAreaSpec & area_spec,
	  const bool & writeout_combinations = false) ;

  /// return an estimate for the number of empty jets -- one uses the
  /// AreaBase one rather than the ActiveArea one (which for which we
  /// do not have the information).
  virtual double n_empty_jets(const Selector & selector) const {
    return ClusterSequenceAreaBase::n_empty_jets(selector);
  }

protected:
  /// does the initialisation and running specific to the passive
  /// areas class
  void _initialise_and_run_1GPA (const JetDefinition & jet_def_in,
                               const GhostedAreaSpec & area_spec,
                               const bool & writeout_combinations = false);

private:

  void _run_1GPA(const GhostedAreaSpec & area_spec);
};




template<class L> ClusterSequence1GhostPassiveArea::ClusterSequence1GhostPassiveArea 
(const std::vector<L> & pseudojets, 
 const JetDefinition & jet_def_in,
 const GhostedAreaSpec & area_spec,
 const bool & writeout_combinations) {

  // transfer the initial jets (type L) into our own array
  _transfer_input_jets(pseudojets);

  // run the clustering for passive areas
  _initialise_and_run_1GPA(jet_def_in, area_spec, writeout_combinations);

}


  
FASTJET_END_NAMESPACE

#endif // __FASTJET_CLUSTERSEQUENCE1GHOSTPASSIVEAREA_HH__
