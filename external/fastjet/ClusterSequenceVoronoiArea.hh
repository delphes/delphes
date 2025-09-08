//FJSTARTHEADER
// $Id$
//
// Copyright (c) 2005-2025, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#ifndef __FASTJET_CLUSTERSEQUENCEVORONOIAREA_HH__
#define __FASTJET_CLUSTERSEQUENCEVORONOIAREA_HH__

#include "fastjet/PseudoJet.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include <memory>
#include <vector>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// @ingroup sec_area_classes
/// \class ClusterSequenceVoronoiArea
/// Like ClusterSequence with computation of the Voronoi jet area
///
/// Handle the computation of Voronoi jet area.
///
/// This class should not be used directly. Rather use
/// ClusterSequenceArea with the appropriate AreaDefinition
class ClusterSequenceVoronoiArea : public ClusterSequenceAreaBase {
public:
  /// template ctor
  /// \param pseudojet              list of jets (template type)
  /// \param jet_def                jet definition
  /// \param effective_Rfact        effective radius
  /// \param writeout_combinations  ??????
  template<class L> ClusterSequenceVoronoiArea
         (const std::vector<L> & pseudojets, 
	  const JetDefinition & jet_def,
	  const VoronoiAreaSpec & spec = VoronoiAreaSpec(),
	  const bool & writeout_combinations = false);
  
  /// default dtor
  ~ClusterSequenceVoronoiArea();

  /// return the area associated with the given jet
  virtual inline double area(const PseudoJet & jet) const FASTJET_OVERRIDE {
    return _voronoi_area[jet.cluster_hist_index()];}

  /// return a 4-vector area associated with the given jet -- strictly
  /// this is not the exact 4-vector area, but rather an approximation
  /// made of sums of centres of all Voronoi cells in jet, each
  /// contributing with a normalisation equal to the area of the cell
  virtual inline PseudoJet area_4vector(const PseudoJet & jet) const FASTJET_OVERRIDE {
    return _voronoi_area_4vector[jet.cluster_hist_index()];}

  /// return the error of the area associated with the given jet
  /// (0 by definition for a voronoi area)
  virtual inline double area_error(const PseudoJet & /*jet*/) const FASTJET_OVERRIDE {
    return 0.0;}

  /// passive area calculator -- to be defined in the .cc file (it will do
  /// the true hard work)
  class VoronoiAreaCalc; 
  

private:  
  /// initialisation of the Voronoi Area
  void _initializeVA();

  std::vector<double> _voronoi_area;  ///< vector containing the result
  std::vector<PseudoJet> _voronoi_area_4vector; ///< vector containing approx 4-vect areas
  VoronoiAreaCalc *_pa_calc;          ///< area calculator
  double _effective_Rfact;            ///< effective radius
};




/// template constructor need to be specified in the header!
//----------------------------------------------------------------------
template<class L> ClusterSequenceVoronoiArea::ClusterSequenceVoronoiArea
(const std::vector<L> &pseudojets, 
 const JetDefinition &jet_def_in,
 const VoronoiAreaSpec & spec,
 const bool & writeout_combinations) :
  _effective_Rfact(spec.effective_Rfact()) {

  // transfer the initial jets (type L) into our own array
  _transfer_input_jets(pseudojets);

  // run the clustering
  _initialise_and_run(jet_def_in,writeout_combinations);

  // the jet clustering's already been done, now worry about areas...
  _initializeVA();
}

FASTJET_END_NAMESPACE

#endif // __FASTJET_CLUSTERSEQUENCEVORONOIAREA_HH__
