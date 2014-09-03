//FJSTARTHEADER
// $Id$
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



#include "fastjet/RectangularGrid.hh"
#include <sstream>
using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

  /// dummy ctor (will give an unusable grid)
RectangularGrid::RectangularGrid()
    : _ymax(-1.0), _ymin(1.0), _requested_drap(-1.0), _requested_dphi(-1.0) {
  // put in nonsensical values for the other variables too, to keep coverity happy
  _ntotal = -1;
  _ngood  = -1;
  _dy = 0.0;
  _dphi = 0.0;
  _cell_area = 0.0;
  _inverse_dy = 0;
  _inverse_dphi = 0;
  _ny   = 0;
  _nphi = 0;
}


int RectangularGrid::tile_index(const PseudoJet & p) const {
  // the code below has seem some degree of optimization: don't change
  // it without testing the speed again

  // new version as of 2014-08-04
  double y_minus_ymin = p.rap() - _ymin;
  if (y_minus_ymin < 0) return -1;
  int iy = int(y_minus_ymin * _inverse_dy); // guaranteed positive, so int is safe
  if (iy >= _ny) return -1;

  // old version: gives a SoftKiller that's about 10% slower on Gavin's retina mac.
  // (though having it hard coded inside SoftKiller returns that advantage)
  // BUT: some comments said that this was a factor of two faster than
  //      something similar to the version above. What is going on?
  // int iy = int(floor( (p.rap() - _ymin) * _inverse_dy ));
  // if (iy < 0 || iy >= _ny) return -1;

  int iphi = int( p.phi() * _inverse_dphi );
  if (iphi == _nphi) iphi = 0; // just in case of rounding errors

  return iy*_nphi + iphi;
}


//----------------------------------------------------------------------
void RectangularGrid::_setup_grid() {
  // initial sanity checks
  assert(_ymax > _ymin);
  assert(_requested_drap > 0);
  assert(_requested_dphi > 0);

  double ny_double = (_ymax-_ymin) / _requested_drap;
  _ny = max(int(ny_double+0.5),1);
  _dy = (_ymax-_ymin) / _ny;
  _inverse_dy = _ny/(_ymax-_ymin);
  
  _nphi = int (twopi / _requested_dphi + 0.5);
  _dphi = twopi / _nphi;
  _inverse_dphi = _nphi/twopi;

  // some sanity checking (could throw a fastjet::Error)
  assert(_ny >= 1 && _nphi >= 1);

  _ntotal = _nphi * _ny;
  //_max_pt.resize(_ntotal);
  _cell_area = _dy * _dphi;

  // if we have a selector, establish which tiles are good;
  // apply the selector to a 4-vector at the tile's centre
  if (_tile_selector.worker()) {
    _is_good.resize(n_tiles());
    _ngood = 0;
    for (int i = 0; i < n_tiles(); i++) {
      int iphi = i % _nphi;
      int irap = i / _nphi;
      double phi = (iphi + 0.5)*_dphi;
      double rap = (irap + 0.5)*_dy + _ymin;
      _is_good[i] = _tile_selector.pass(PtYPhiM(1.0, rap, phi));
      if (_is_good[i]) _ngood++;
    }
  } else {
    _ngood = n_tiles();
  }
}

//----------------------------------------------------------------------
string RectangularGrid::description() const {
  if (! is_initialised())
    return "Uninitialised rectangular grid";

  ostringstream oss;
  oss << "rectangular grid with rapidity extent " << _ymin << " < rap < " << _ymax
      << ", tile size drap x dphi = " << _dy << " x " << _dphi;

  if (_tile_selector.worker()) {
    oss << ", good tiles are those that pass selector " <<  _tile_selector.description();
  }
  return oss.str();
}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
