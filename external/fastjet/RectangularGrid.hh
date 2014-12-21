#ifndef __FASTJET_RECTANGULARGRID_HH__
#define __FASTJET_RECTANGULARGRID_HH__

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

#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


//----------------------------------------------------------------------
/// Class to indicate generic structure of tilings
class TilingBase {
public:
  virtual ~TilingBase() {}

  /// returns the index of the tile in which p is located, or -1 if p
  /// is outside the tiling region
  virtual int tile_index(const PseudoJet & p) const = 0;

  /// returns the total number of tiles in the tiling; valid tile
  /// indices run from 0 ... n_tiles()-1;
  virtual int n_tiles() const = 0;

  /// returns the number of tiles that are "good"; i.e. there is scope
  /// for having tiles that, for whatever reason, should be ignored;
  /// there are situations in which having "non-good" tiles may be the
  /// simplest mechanism to obtain a tiling with holes in it
  virtual int n_good_tiles() const {return n_tiles();}

  /// returns whether a given tile is good
  virtual bool tile_is_good(int /* itile */) const {return true;}

  /// returns whether all tiles are good
  virtual bool all_tiles_good() const {return n_good_tiles() == n_tiles();}

  /// returns true if all tiles have the same area
  virtual bool all_tiles_equal_area() const {return true;}

  /// returns the area of tile itile. Here with a default
  /// implementation to return mean_tile_area(), consistent with the
  /// fact that all_tiles_equal_area() returns true.
  virtual double tile_area(int /* itile */) const {return mean_tile_area();}

  /// returns the mean area of the tiles.
  virtual double mean_tile_area() const = 0;

  /// returns a string to describe the tiling
  virtual std::string description() const = 0;

  /// returns true if the Tiling structure is in a suitably initialised state
  virtual bool is_initialised() const = 0;
  bool is_initialized() const {return is_initialised();}
};

//----------------------------------------------------------------------
/// Class that holds a generic rectangular tiling
class RectangularGrid : public TilingBase {
public:
  /// ctor with simple initialisation
  ///  \param rapmax     the maximal absolute rapidity extent of the grid
  ///  \param cell_size  the grid spacing (equivalently, cell size)
  RectangularGrid(double rapmax_in, double cell_size) :
      _ymax(rapmax_in), _ymin(-rapmax_in), 
      _requested_drap(cell_size), _requested_dphi(cell_size) {
    _setup_grid();
  }

  /// ctor with more control over initialisation
  ///  \param rapmin         the minimum rapidity extent of the grid
  ///  \param rapmax         the maximum rapidity extent of the grid
  ///  \param drap           the grid spacing in rapidity
  ///  \param dphi           the grid spacing in azimuth
  ///  \param tile_selector  optional (geometric) selector to specify 
  ///                        which tiles are good; a tile is good if
  ///                        a massless 4-vector at the center of the tile passes
  ///                        the selection
  RectangularGrid(double rapmin_in, double rapmax_in, double drap_in, double dphi_in,
                  Selector tile_selector = Selector()) 
    : _ymax(rapmax_in), _ymin(rapmin_in), 
      _requested_drap(drap_in), _requested_dphi(dphi_in),
      _tile_selector(tile_selector)
  {
    _setup_grid();
  }

  /// dummy ctor (will give an unusable grid)
  RectangularGrid();

  virtual int n_tiles() const {return _ntotal;}

  virtual int n_good_tiles() const {return _ngood;}

  // this was being kept inline, but it seems to make little
  // difference whether it is or not (at least on Gavin's mac)
  virtual int tile_index(const PseudoJet & p) const;

  /// returns whether a given tile is good
  // tested in "issue" 2014-08-08-testing-rect-grid
  virtual bool tile_is_good(int itile) const {return _tile_selector.worker() ? _is_good[itile] : true;}

  /// returns the area of tile itile.
  virtual double tile_area(int /* itile */) const {return mean_tile_area();}

  /// returns the mean area of tiles.
  virtual double mean_tile_area() const {return _dphi*_dy;};

  /// returns a textual description of the grid
  virtual std::string description() const;
  
  /// returns the minimum rapidity extent of the grid
  double rapmin() const {return _ymin;}
  /// returns the maxmium rapidity extent of the grid
  double rapmax() const {return _ymax;}
  /// returns the spacing of the grid in rapidity
  double drap()   const {return _dy;}
  /// returns the spacing of the grid in azimuth
  double dphi()   const {return _dphi;}

  /// returns true if the grid is in a suitably initialised state
  virtual bool is_initialised() const {return _ntotal > 0;}

private:
  void _setup_grid();
  
  // information about the requested grid
  double _ymax, _ymin;  ///< maximal and minimal rapidity coverage of the grid
  double _requested_drap; ///< requested rapidity spacing
  double _requested_dphi; ///< requested phi spacing

  // information about the actual grid
  double _dy, _dphi, _cell_area, _inverse_dy, _inverse_dphi;
  int _ny, _nphi, _ntotal;
  int _ngood;

  // a tile selector
  Selector _tile_selector;
  // a cached 
  std::vector<bool> _is_good;
  
};

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __FASTJET_RECTANGULARGRID_HH__
