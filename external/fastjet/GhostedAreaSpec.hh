//FJSTARTHEADER
// $Id: GhostedAreaSpec.hh 3433 2014-07-23 08:17:03Z salam $
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


#ifndef __FASTJET_GHOSTEDAREASPEC_HH__
#define __FASTJET_GHOSTEDAREASPEC_HH__

#include<vector>
#include<string>
#include "fastjet/PseudoJet.hh"
#include "fastjet/internal/BasicRandom.hh"
#include "fastjet/Selector.hh"
#include "fastjet/LimitedWarning.hh"

// 
#define STATIC_GENERATOR 1

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// namespace to hold default parameters for the active area spec
namespace gas {
  const double def_ghost_maxrap  = 6.0;
  const int    def_repeat        = 1;
  const double def_ghost_area    = 0.01;
  const double def_grid_scatter  = 1.0;
  const double def_pt_scatter    = 0.1;
  const double def_mean_ghost_pt = 1e-100;
}

//----------------------------------------------------------------------
/// @ingroup area_classes
/// \class GhostedAreaSpec
/// Parameters to configure the computation of jet areas using ghosts
///
/// Class that defines the parameters that go into the measurement
/// of active jet areas.
class GhostedAreaSpec {
public:
  /// default constructor
  GhostedAreaSpec(): _ghost_maxrap (gas::def_ghost_maxrap), 
		     _ghost_rap_offset(0.0),
                    _repeat       (gas::def_repeat), 
                    _ghost_area   (gas::def_ghost_area), 
                    _grid_scatter (gas::def_grid_scatter), 
                    _pt_scatter   (gas::def_pt_scatter), 
                    _mean_ghost_pt(gas::def_mean_ghost_pt),
                    _fj2_placement(false) {_initialize();}
  
  /// explicit constructor
  explicit GhostedAreaSpec(double ghost_maxrap_in, 
                          int    repeat_in        = gas::def_repeat,
                          double ghost_area_in    = gas::def_ghost_area,   
                          double grid_scatter_in  = gas::def_grid_scatter, 
                          double pt_scatter_in    = gas::def_pt_scatter,   
                          double mean_ghost_pt_in = gas::def_mean_ghost_pt
                          ): 
    _ghost_maxrap(ghost_maxrap_in), 
    _ghost_rap_offset(0.0),
    _repeat(repeat_in), 
    _ghost_area(ghost_area_in), 
    _grid_scatter(grid_scatter_in),  
    _pt_scatter(pt_scatter_in), 
    _mean_ghost_pt(mean_ghost_pt_in),
    _fj2_placement(false) {_initialize();}

  /// explicit constructor
  explicit GhostedAreaSpec(double ghost_minrap_in, 
			   double ghost_maxrap_in, 
                           int    repeat_in        = gas::def_repeat,
                           double ghost_area_in    = gas::def_ghost_area,   
                           double grid_scatter_in  = gas::def_grid_scatter, 
                           double pt_scatter_in    = gas::def_pt_scatter,   
                           double mean_ghost_pt_in = gas::def_mean_ghost_pt
                          ): 
    _ghost_maxrap    (0.5*(ghost_maxrap_in - ghost_minrap_in)), 
    _ghost_rap_offset(0.5*(ghost_maxrap_in + ghost_minrap_in)),
    _repeat(repeat_in), 
    _ghost_area(ghost_area_in), 
    _grid_scatter(grid_scatter_in),  
    _pt_scatter(pt_scatter_in), 
    _mean_ghost_pt(mean_ghost_pt_in),
    _fj2_placement(false) {_initialize();}


  /// constructor based on a Selector
  explicit GhostedAreaSpec(const Selector & selector,
                           int    repeat_in        = gas::def_repeat,
                           double ghost_area_in    = gas::def_ghost_area,   
                           double grid_scatter_in  = gas::def_grid_scatter, 
                           double pt_scatter_in    = gas::def_pt_scatter,   
                           double mean_ghost_pt_in = gas::def_mean_ghost_pt
			   );


  /// does the initialization of actual ghost parameters
  void _initialize();

  // for accessing values set by the user
  inline double ghost_rapmax () const {return _ghost_maxrap;}
  inline double ghost_maxrap () const {return _ghost_maxrap;}
  inline double ghost_etamax () const {return _ghost_maxrap;}
  inline double ghost_maxeta () const {return _ghost_maxrap;}
  inline double ghost_area   () const {return _ghost_area   ;}
  inline double grid_scatter () const {return _grid_scatter;}
  inline double pt_scatter   () const {return _pt_scatter  ;}
  inline double mean_ghost_pt() const {return _mean_ghost_pt  ;}
  inline int    repeat       () const {return _repeat      ;}
  inline bool   fj2_placement() const {return _fj2_placement;}

  inline double kt_scatter   () const {return _pt_scatter  ;}
  inline double mean_ghost_kt() const {return _mean_ghost_pt  ;}

  // for accessing values 
  inline double actual_ghost_area() const {return _actual_ghost_area;}
  inline int    n_ghosts()          const {return _n_ghosts;}

  // when explicitly modifying values, sometimes call the initializer
  inline void set_ghost_area   (double val) {_ghost_area    = val; _initialize();}
  inline void set_ghost_rapmax (double val) {_ghost_maxrap = val; _initialize();}
  inline void set_ghost_maxrap (double val) {_ghost_maxrap = val; _initialize();}
  inline void set_ghost_etamax (double val) {_ghost_maxrap = val; _initialize();}
  inline void set_ghost_maxeta (double val) {_ghost_maxrap = val; _initialize();}
  inline void set_grid_scatter (double val) {_grid_scatter   = val; }
  inline void set_pt_scatter   (double val) {_pt_scatter     = val; }
  inline void set_mean_ghost_pt(double val) {_mean_ghost_pt  = val; }
  inline void set_repeat       (int    val) {_repeat         = val; }

  inline void set_kt_scatter   (double val) {_pt_scatter     = val; }
  inline void set_mean_ghost_kt(double val) {_mean_ghost_pt  = val; }

  /// if val is true, set ghost placement as it was in FastJet 2.X. The
  /// main differences between FJ2 and FJ3 ghost placement are
  ///
  ///  - in FJ2 the rapidity spacing was
  ///    ceil((maxrap-minrap)/sqrt(area)), while in FJ3 it is
  ///    int((maxrap-minrap)/sqrt(area) + 0.5) [similarly for phi].
  ///    The FJ3 option offers more stability when trying to specify a
  ///    spacing that exactly fits the extent.
  ///
  /// - in FJ2, the ghosts are placed at the corners of grid cells
  ///   (i.e. extending up to maxrap), while in FJ3 they are placed at
  ///   the centres of grid cells (i.e. extending roughly up to
  ///   maxrap-sqrt(area)). The FJ2 behaviour effectively skews the
  ///   total area coverage when maxrap is small, by an amount
  ///   sqrt(area)/(2*maxrap).
  ///
  /// FJ2 placement is now deprecated.
  void set_fj2_placement(bool  val);

  /// return nphi (ghosts layed out (-nrap, 0..nphi-1), (-nrap+1,0..nphi-1),
  /// ... (nrap,0..nphi-1)
  inline int nphi() const {return _nphi;}
  inline int nrap() const {return _nrap;}

  /// get all relevant information about the status of the 
  /// random number generator, so that it can be reset subsequently
  /// with set_random_status.
  inline void get_random_status(std::vector<int> & __iseed) const {
    _random_generator.get_status(__iseed);}

  /// set the status of the random number generator, as obtained
  /// previously with get_random_status. Note that the random
  /// generator is a static member of the class, i.e. common to all
  /// instances of the class --- so if you modify the random for this
  /// instance, you modify it for all instances.
  inline void set_random_status(const std::vector<int> & __iseed) {
    _random_generator.set_status(__iseed);}
  
  inline void checkpoint_random() {get_random_status(_random_checkpoint);}
  inline void restore_checkpoint_random() {set_random_status(_random_checkpoint);}

  /// for a summary
  std::string description() const;

  /// push a set of ghost 4-momenta onto the back of the vector of
  /// PseudoJets
  void add_ghosts(std::vector<PseudoJet> & ) const;

  /// very deprecated public access to a random number 
  /// from the internal generator
  inline double random_at_own_risk() const {return _our_rand();}
  /// very deprecated public access to the generator itself
  inline BasicRandom<double> & generator_at_own_risk() const {
    return _random_generator;}

private:
  
  // quantities that determine nature and distribution of ghosts
  double _ghost_maxrap;
  double _ghost_rap_offset;
  int    _repeat      ;
  double _ghost_area   ;  
  double _grid_scatter;
  double _pt_scatter  ;
  double _mean_ghost_pt;
  bool   _fj2_placement;

  Selector _selector;

  // derived quantities
  double _actual_ghost_area, _dphi, _drap;
  int    _n_ghosts, _nphi, _nrap;


  std::vector<int> _random_checkpoint;
  static BasicRandom<double> _random_generator;
  //mutable BasicRandom<double> _random_generator;

  static LimitedWarning _warn_fj2_placement_deprecated;

  inline double _our_rand() const {return _random_generator();}
  
};

/// just provide a typedef for backwards compatibility with programs
/// based on versions 2.0 and 2.1 of fastjet. Since there is no
/// easy way of telling people this is deprecated at compile or run
/// time, we should be careful before removing this in the future.
typedef GhostedAreaSpec ActiveAreaSpec;


FASTJET_END_NAMESPACE

#endif // __FASTJET_GHOSTEDAREASPEC_HH__
