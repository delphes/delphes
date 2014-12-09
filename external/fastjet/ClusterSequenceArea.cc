//FJSTARTHEADER
// $Id: ClusterSequenceArea.cc 3433 2014-07-23 08:17:03Z salam $
//
// Copyright (c) 2006-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#include "fastjet/ClusterSequenceArea.hh"

FASTJET_BEGIN_NAMESPACE

LimitedWarning ClusterSequenceArea::_range_warnings;
LimitedWarning ClusterSequenceArea::_explicit_ghosts_repeats_warnings;

/// print a warning if the range is unsuitable for the current
/// calculation of the area (e.g. because ghosts do not extend
/// far enough).
void ClusterSequenceArea::_warn_if_range_unsuitable(const Selector & selector) const {
  _check_selector_good_for_median(selector);

  bool no_ghosts = (_area_def.area_type() == voronoi_area)
    || (_area_def.area_type() == passive_area
        && jet_def().jet_algorithm() == kt_algorithm);
  if (! no_ghosts) {
    double rapmin, rapmax;
    selector.get_rapidity_extent(rapmin, rapmax);
    if (rapmin < -_area_def.ghost_spec().ghost_maxrap()+0.95*jet_def().R() ||
        rapmax >  _area_def.ghost_spec().ghost_maxrap()-0.95*jet_def().R()) {
      _range_warnings.warn("rapidity range for median (rho) extends beyond +-(ghost_maxrap - 0.95*R); this is likely to cause the results to be unreliable; safest option is to increase ghost_maxrap in the area definition");
    }
  }
}


FASTJET_END_NAMESPACE
