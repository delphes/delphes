//STARTHEADER
// $Id: CDFMidPointPlugin.hh 2758 2011-11-24 08:31:58Z soyez $
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

#ifndef __CDFMIDPOINTPLUGIN_HH__
#define __CDFMIDPOINTPLUGIN_HH__

#include "fastjet/JetDefinition.hh"

// questionable whether this should be in fastjet namespace or not...

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
//
/// @ingroup plugins
/// \class CDFMidPointPlugin
/// Implementation of the MidPoint algorithm from CDF (plugin for
/// fastjet-v2.1 upwards)
///
/// A plugin for fastjet-v2.1 that provides an interface to the CDF
/// midpoint algorithm
///
/// CDFMidPointPlugin is a plugin for fastjet (v2.1 upwards) that
/// provides an interface to the CDF version of Run-II iterative cone
/// algorithm with midpoint seeds (also known as the Iterative Legacy
/// Cone Algorithm, ILCA).
///
/// The CDF code has been taken from Joey Huston's webpage
/// http://www.pa.msu.edu/~huston/Les_Houches_2005/Les_Houches_SM.html
///
/// Note that the CDF midpoint code contains options that go beyond
/// those described in the Tevatron run-II document (hep-ex/0005012),
/// notably search-cones, as described in hep-ph/0111434, and
/// midpoints bewteen multiplets of stable cones.
///
/// Additionally, the version of the CDF midpoint code distributed
/// here has been modified by the FastJet authors, so as to allow one
/// to choose the scale used in the split-merge step.
//
//----------------------------------------------------------------------
class CDFMidPointPlugin : public JetDefinition::Plugin {
public:
  /// the choice of scale to be used in the split-merge step
  // NB: just replicates what we've added to the CDF midpoint code
  enum SplitMergeScale {SM_pt, SM_Et, SM_mt, SM_pttilde};

  ///
  /// A CDFMidPointPlugin constructor that looks like the one provided
  /// by CDF. Its arguments should have the following meaning:
  ///
  /// - seed_threshold: minimum pt for a particle to be considered 
  ///   a seed of the iteration.
  ///
  /// - cone_radius: standard meaning
  ///
  /// - cone_area_fraction: stable-cones are searched for with a
  ///   radius Rsearch = R * sqrt(cone_area_fraction), and then
  ///   expanded to size R afterwards; note (hep-ph/0610012) that this
  ///   introduces IR unsafety at NLO for X+2-jet observables (where X
  ///   any hard object).
  ///
  /// - max_pair_size: "midpoints" can be added between pairs of
  ///   stable cones, triplets of stable cones, etc.; max_pair_size
  ///   indicates the maximum number of stable cones that are
  ///   assembled when adding midpoints.
  ///
  /// - max_iterations: the maximum number of iterations to carry out
  ///   when looking for a stable cone.
  ///
  /// - overlap_threshold: if
  ///     (overlapping_Et)/(Et_of_softer_protojet) < overlap_threshold,
  ///   overlapping jets are split, otherwise they are merged.
  ///
  /// - sm_scale: a choice for the scale to be used in the split-merge
  ///   step (both for ordering the momenta and quantifying the
  ///   overlap); the three options are
  ///
  ///    . SM_pt: pt (default -- source of small IR safety issue in purely
  ///      hadronic events)
  ///
  ///    . SM_Et: Et (not boost invariant, reduces to mt at zero rapidity and
  ///      to pt and infinite rapidity)
  ///
  ///    . SM_mt: transverse mass = sqrt(m^2+pt^2)
  ///
  CDFMidPointPlugin (
                     double seed_threshold_in     ,	 
		     double cone_radius_in        ,
		     double cone_area_fraction_in ,
		     int    max_pair_size_in      ,
		     int    max_iterations_in     ,
		     double overlap_threshold_in  ,
                     SplitMergeScale sm_scale_in = SM_pt) :
    _seed_threshold     (seed_threshold_in     ),    
    _cone_radius        (cone_radius_in        ),
    _cone_area_fraction (cone_area_fraction_in ),
    _max_pair_size      (max_pair_size_in      ),
    _max_iterations     (max_iterations_in     ),
    _overlap_threshold  (overlap_threshold_in  ),
    _sm_scale           (sm_scale_in)             {}

  /// a compact constructor
  ///
  /// NB: as of version 2.4, the default value for the
  /// overlap_threshold threshold has been removed, to avoid
  /// misleading people into using the value of 0.5 without thinking,
  /// which is known to have adverse effects in high-noise
  /// environments. A recommended value is 0.75.
  CDFMidPointPlugin (double   cone_radius_in, 
		     double   overlap_threshold_in,// = 0.5, 
		     double   seed_threshold_in = 1.0,	     
		     double   cone_area_fraction_in = 1.0) : 
    _seed_threshold     (seed_threshold_in     ),    
    _cone_radius        (cone_radius_in        ),
    _cone_area_fraction (cone_area_fraction_in ),
    _max_pair_size      (2                     ),
    _max_iterations     (100                   ),
    _overlap_threshold  (overlap_threshold_in  ),
    _sm_scale           (SM_pt)                {}


  // some functions to return info about parameters
  double seed_threshold     () const {return _seed_threshold     ;}
  double cone_radius        () const {return _cone_radius        ;}
  double cone_area_fraction () const {return _cone_area_fraction ;}
  int    max_pair_size      () const {return _max_pair_size      ;}
  int    max_iterations     () const {return _max_iterations     ;}
  double overlap_threshold  () const {return _overlap_threshold  ;}


  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;
  /// the plugin mechanism's standard way of accessing the jet radius
  virtual double R() const {return cone_radius();}
                      
private:

  double _seed_threshold    ;
  double _cone_radius       ;
  double _cone_area_fraction;
  int    _max_pair_size     ;
  int    _max_iterations    ;
  double _overlap_threshold ;
  SplitMergeScale _sm_scale ;

  static bool _first_time;

  /// print a banner for reference to the 3rd-party code
  void _print_banner(std::ostream *ostr) const;
};

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif // __CDFMIDPOINTPLUGIN_HH__
