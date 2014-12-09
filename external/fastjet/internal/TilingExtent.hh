#ifndef __FASTJET_TILINGEXTENT_HH__
#define __FASTJET_TILINGEXTENT_HH__

#include "fastjet/ClusterSequence.hh"

//FJSTARTHEADER
// $Id: TilingExtent.hh 3433 2014-07-23 08:17:03Z salam $
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

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
/// class to perform a fast analysis of the appropriate rapidity range
/// in which to perform tiling
class TilingExtent {
public:
  /// constructor that takes a ClusterSequence in a state where the
  /// initial particles have been set up, but before clustering has
  /// started.
  TilingExtent(ClusterSequence & cs);

  /// returns the suggested minimum rapidity for the tiling
  double minrap() const {return _minrap;}
  /// returns the suggested maximum rapidity for the tiling
  double maxrap() const {return _maxrap;}

  /// internally, the class bins the particle multiplicity versus
  /// rapidity, in bins of size 1 running roughly from minrap to maxrap
  /// (including overflows); this function returns the sum of squares
  /// of bin contents, which may be informative for deciding strategy
  /// choices.
  double sum_of_binned_squared_multiplicity() const {return _cumul2;}

private:
  double _minrap, _maxrap, _cumul2;

  /// attempts to calculate a sensible rapidity extent for the tiling
  void _determine_rapidity_extent(const std::vector<PseudoJet> & particles);
};

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif // __FASTJET_TILINGEXTENT_HH__
