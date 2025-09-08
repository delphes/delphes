#ifndef __GRID_MEDIAN_BACKGROUND_ESTIMATOR_HH__
#define __GRID_MEDIAN_BACKGROUND_ESTIMATOR_HH__

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


#include "fastjet/tools/BackgroundEstimatorBase.hh"
#include "fastjet/RectangularGrid.hh"


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// @ingroup tools_background
/// \class GridMedianBackgroundEstimator
/// 
/// Background Estimator based on the median pt/area of a set of grid
/// cells. 
///
/// Description of the method:
///   This background estimator works by projecting the event onto a
///   grid in rapidity and azimuth. In each grid cell, the scalar pt
///   sum of the particles in the cell is computed. The background
///   density is then estimated by the median of (scalar pt sum/cell
///   area) for all cells.
///
/// Parameters:
///   The class takes 2 arguments: the absolute rapidity extent of the 
///   cells and the size of the grid cells. Note that the size of the cell
///   will be adjusted in azimuth to satisfy the 2pi periodicity and
///   in rapidity to match the requested rapidity extent.
///
/// Rescaling:
///   It is possible to use a rescaling profile. In this case, the
///   profile needs to be set before setting the particles and it will
///   be applied to each particle (i.e. not to each cell). 
///   Note also that in this case one needs to call rho(jet) instead of
///   rho() [Without rescaling, they are identical]
///
class GridMedianBackgroundEstimator : public BackgroundEstimatorBase
                                    , public RectangularGrid
{

public:
  /// @name  constructors and destructors
  //\{
  //----------------------------------------------------------------
  ///   \param ymax   maximal absolute rapidity extent of the grid
  ///   \param requested_grid_spacing   size of the grid cell. The
  ///            "real" cell size could differ due e.g. to the 2pi
  ///             periodicity in azimuthal angle (size, not area)
  GridMedianBackgroundEstimator(double ymax, double requested_grid_spacing) :
    RectangularGrid(ymax, requested_grid_spacing),
    _enable_rho_m(true) {} 

  //----------------------------------------------------------------
  /// Constructor based on a user's fully specified RectangularGrid
  GridMedianBackgroundEstimator(const RectangularGrid & grid) :
    RectangularGrid(grid), _enable_rho_m(true) {
    if (!RectangularGrid::is_initialised()) 
      throw Error("attempt to construct GridMedianBackgroundEstimator with uninitialised RectangularGrid");
  }    

  //---------------------------------------------------------------- 
  /// Constructor with the explicit parameters for the underlying
  /// RectangularGrid
  ///
  ///  \param rapmin         the minimum rapidity extent of the grid
  ///  \param rapmax         the maximum rapidity extent of the grid
  ///  \param drap           the grid spacing in rapidity
  ///  \param dphi           the grid spacing in azimuth
  ///  \param tile_selector  optional (geometric) selector to specify 
  ///                        which tiles are good; a tile is good if
  ///                        a massless 4-vector at the center of the tile passes
  ///                        the selection
  GridMedianBackgroundEstimator(double rapmin_in, double rapmax_in, double drap_in, double dphi_in,
                                Selector tile_selector = Selector()) :
    RectangularGrid(rapmin_in, rapmax_in, drap_in, dphi_in, tile_selector), _enable_rho_m(true) {}

  //\}


  /// @name setting a new event
  //\{
  //----------------------------------------------------------------

  /// tell the background estimator that it has a new event, composed
  /// of the specified particles.
  void set_particles(const std::vector<PseudoJet> & particles) FASTJET_OVERRIDE;

  /// determine whether the automatic calculation of rho_m and sigma_m
  /// is enabled (by default true)
  void set_compute_rho_m(bool enable){ _enable_rho_m = enable; }

  //\}

  /// return a pointer to a copy of this BGE; the user is responsible
  /// for eventually deleting the resulting object.
  BackgroundEstimatorBase * copy() const FASTJET_OVERRIDE {
    return new GridMedianBackgroundEstimator(*this);
  };



  /// @name  retrieving fundamental information
  //\{
  //----------------------------------------------------------------
  /// get the full set of background properties
  BackgroundEstimate estimate() const FASTJET_OVERRIDE;
  
  /// get the full set of background properties for a given reference jet
  BackgroundEstimate estimate(const PseudoJet &jet) const FASTJET_OVERRIDE;

  /// returns rho, the median background density per unit area
  double rho() const FASTJET_OVERRIDE;

  /// returns sigma, the background fluctuations per unit area; must be
  /// multipled by sqrt(area) to get fluctuations for a region of a
  /// given area.
  double sigma() const FASTJET_OVERRIDE;

  /// returns rho, the background density per unit area, locally at the
  /// position of a given jet. Note that this is not const, because a
  /// user may then wish to query other aspects of the background that
  /// could depend on the position of the jet last used for a rho(jet)
  /// determination.
  double rho(const PseudoJet & jet) FASTJET_OVERRIDE;

  /// returns sigma, the background fluctuations per unit area, locally at
  /// the position of a given jet. As for rho(jet), it is non-const.
  double sigma(const PseudoJet & jet) FASTJET_OVERRIDE;

  /// returns true if this background estimator has support for
  /// determination of sigma
  bool has_sigma() const FASTJET_OVERRIDE {return true;}

  //-----------------------------------------------------------------
  /// Returns rho_m, the purely longitudinal, particle-mass-induced
  /// component of the background density per unit area
  double rho_m() const FASTJET_OVERRIDE;

  /// returns sigma_m, a measure of the fluctuations in the purely
  /// longitudinal, particle-mass-induced component of the background
  /// density per unit area; must be multipled by sqrt(area) to get
  /// fluctuations for a region of a given area.
  double sigma_m() const FASTJET_OVERRIDE;

  /// Returns rho_m locally at the jet position. As for rho(jet), it is non-const.
  double rho_m(const PseudoJet & jet) FASTJET_OVERRIDE;

  /// Returns sigma_m locally at the jet position. As for rho(jet), it is non-const.
  double sigma_m(const PseudoJet & jet) FASTJET_OVERRIDE;

  /// Returns true if this background estimator has support for
  /// determination of rho_m.
  ///
  /// Note that support for sigma_m is automatic if one has sigma and
  /// rho_m support.
  bool has_rho_m() const FASTJET_OVERRIDE {return _enable_rho_m;}


  /// returns the area of the grid cells (all identical, but
  /// referred to as "mean" area for uniformity with JetMedianBGE).
  double mean_area() const {return mean_tile_area();}
  //\}

  /// @name configuring the behaviour
  //\{
  //----------------------------------------------------------------

  /// Set a pointer to a class that calculates the rescaling factor as
  /// a function of the jet (position). Note that the rescaling factor
  /// is used both in the determination of the "global" rho (the pt/A
  /// of each jet is divided by this factor) and when asking for a
  /// local rho (the result is multiplied by this factor).
  ///
  /// The BackgroundRescalingYPolynomial class can be used to get a
  /// rescaling that depends just on rapidity.
  ///
  /// Note that this has to be called BEFORE any attempt to do an
  /// actual computation
  ///
  /// The same profile will be used for both pt and mt (this is
  /// probabaly a good approximation since the particle density
  /// changes is what dominates the rapidity profile)
  virtual void set_rescaling_class(const FunctionOfPseudoJet<double> * rescaling_class) FASTJET_OVERRIDE;

  //\}

  /// @name description
  //\{
  //----------------------------------------------------------------

  /// returns a textual description of the background estimator
  std::string description() const FASTJET_OVERRIDE;

  //\}


private:

  /// verify that particles have been set and throw an error if not
  void verify_particles_set() const;

  // information about the event
  //std::vector<double> _scalar_pt;
  //double _rho, _sigma, _rho_m, _sigma_m;
  BackgroundEstimate _cached_estimate;
  //bool _has_particles;
  bool _enable_rho_m;

  // various warnings to inform people of potential dangers
  LimitedWarning _warning_rho_of_jet;
  LimitedWarning _warning_rescaling;
};

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __GRID_MEDIAN_BACKGROUND_ESTIMATOR_HH__
