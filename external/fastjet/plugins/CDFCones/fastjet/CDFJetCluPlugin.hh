//STARTHEADER
// $Id: CDFJetCluPlugin.hh 2758 2011-11-24 08:31:58Z soyez $
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

#ifndef __CDFJETCLUPLUGIN_HH__
#define __CDFJETCLUPLUGIN_HH__

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include <map>

// questionable whether this should be in fastjet namespace or not...

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// @ingroup plugins
/// \class CDFJetCluPlugin
/// Implementation of the JetClu algorithm from CDF (plugin for
/// fastjet-v2.1 upwards)
class CDFJetCluPlugin : public JetDefinition::Plugin {
public:
  /// a compact constructor
  CDFJetCluPlugin (double   cone_radius_in, 
		   double   overlap_threshold_in, 
		   double   seed_threshold_in = 1.0,
		   int      iratch_in = 1) : 
    _seed_threshold    ( seed_threshold_in    ),    
    _cone_radius       ( cone_radius_in       ),
    _adjacency_cut     (   2                  ),
    _max_iterations    ( 100                  ),
    _iratch            ( iratch_in            ),
    _overlap_threshold ( overlap_threshold_in )  {}

  /// a constructor that looks like the one provided by CDF
  CDFJetCluPlugin (double seed_threshold_in   ,	 
		   double cone_radius_in      ,
		   int    adjacency_cut_in    ,
		   int    max_iterations_in   ,
		   int    iratch_in           ,
		   double overlap_threshold_in) :
    _seed_threshold    (seed_threshold_in    ),    
    _cone_radius       (cone_radius_in       ),
    _adjacency_cut     (adjacency_cut_in     ),
    _max_iterations    (max_iterations_in    ),
    _iratch            (iratch_in            ),
    _overlap_threshold (overlap_threshold_in )  {}

  // some functions to return info about parameters
  double seed_threshold    () const {return _seed_threshold    ;}
  double cone_radius       () const {return _cone_radius       ;}
  int    adjacency_cut     () const {return _adjacency_cut     ;}
  int    max_iterations    () const {return _max_iterations    ;}
  int    iratch            () const {return _iratch            ;}
  double overlap_threshold () const {return _overlap_threshold ;}


  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;
  /// the plugin mechanism's standard way of accessing the jet radius
  virtual double R() const {return cone_radius();}
                      
private:

  double _seed_threshold   ;
  double _cone_radius      ;
  int    _adjacency_cut    ;
  int    _max_iterations   ;
  int    _iratch           ;
  double _overlap_threshold;

  /// given a jet try inserting its energy into the map -- if that
  /// energy entry already exists, modify the jet infinitesimally so
  /// as ensure that the jet energy is unique
  void _insert_unique (PseudoJet & jet, std::map<double,int> & jetmap) const;

  static bool _first_time;

  /// print a banner for reference to the 3rd-party code
  void _print_banner(std::ostream *ostr) const;
};

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif // __CDFJETCLUPLUGIN_HH__
