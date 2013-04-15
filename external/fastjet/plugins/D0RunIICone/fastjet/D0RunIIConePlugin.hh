#ifndef __D0RUNIICONEPLUGIN_HH__
#define __D0RUNIICONEPLUGIN_HH__

//STARTHEADER
// $Id: D0RunIIConePlugin.hh 2761 2011-11-24 13:54:05Z soyez $
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

#include "fastjet/JetDefinition.hh"

// questionable whether this should be in fastjet namespace or not...

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
//
/// @ingroup plugins
/// \class D0RunIIConePlugin
/// Implementation of the D0 Run II Cone (plugin for fastjet v2.1 upwards)
/// 
/// D0RunIIConePlugin is a plugin for fastjet (v2.1 upwards) that
/// provides an interface to the D0 version of Run-II iterative cone
/// algorithm with midpoint seeds (also known as the Iterative Legacy
/// Cone Algorithm, ILCA).
///
/// The D0 code has been taken from Lars Sonnenschein's web-space
/// http://www-d0.fnal.gov/~sonne/D0RunIIcone.tgz
///
/// The version of the D0 Run II code distributed
/// here has been modified by the FastJet authors, so as to provide
/// access to the contents of the jets (as is necessary for the
/// plugin). This does not modify the results of the clustering.
//
//----------------------------------------------------------------------
class D0RunIIConePlugin : public JetDefinition::Plugin {
public:

  //
  /// A D0RunIIConePlugin constructor which sets the "free" parameters of the
  /// algorithm:
  ///
  ///  - the cone_radius has the usual meaning
  ///
  ///  - the min_jet_Et causes cones to be discarded at if at any
  ///    iteration they have pt < Et_min_ratio * min_jet_Et. Two
  ///    values have been used by D0 for min_jet_Et: 8 GeV in earlier
  ///    Run II publicatinos, 6 GeV in later publications
  ///
  ///  - split_ratio is equivalent to the overlap threshold during the split/merge step. 
  ///    Default: 0.5.
  ///
  /// The remaining parameters of the algorithm are not to be modified if the algorithm
  /// is to correspond to the one actually used by D0.
  //
  D0RunIIConePlugin (double cone_radius_in, 
                     double min_jet_Et_in , 
                     double split_ratio_in = _DEFAULT_split_ratio) :
    _cone_radius            (cone_radius_in                  ),
    _min_jet_Et             (min_jet_Et_in                   ),
    _split_ratio            (split_ratio_in                  ),
    _far_ratio              (_DEFAULT_far_ratio              ),
    _Et_min_ratio           (_DEFAULT_Et_min_ratio           ),
    _kill_duplicate         (_DEFAULT_kill_duplicate         ),
    _duplicate_dR           (_DEFAULT_duplicate_dR           ),
    _duplicate_dPT          (_DEFAULT_duplicate_dPT          ),
    _search_factor          (_DEFAULT_search_factor          ),
    _pT_min_leading_protojet(_DEFAULT_pT_min_leading_protojet),
    _pT_min_second_protojet (_DEFAULT_pT_min_second_protojet ),
    _merge_max              (_DEFAULT_merge_max              ),
    _pT_min_nomerge         (_DEFAULT_pT_min_nomerge         )
  {
    // nothing to be done here!
  }

  // some functions to return info about parameters
  inline double cone_radius            () const { return _cone_radius            ;} //= 0.5;
  inline double min_jet_Et             () const { return _min_jet_Et             ;} //= 8.0;
  inline double split_ratio            () const { return _split_ratio            ;} //= 0.5;
  inline double far_ratio              () const { return _far_ratio              ;} // =0.5;
  inline double Et_min_ratio           () const { return _Et_min_ratio           ;} // =0.5;
  inline bool   kill_duplicate         () const { return _kill_duplicate         ;} // =true;
  inline double duplicate_dR           () const { return _duplicate_dR           ;} // =0.005; 
  inline double duplicate_dPT          () const { return _duplicate_dPT          ;} // =0.01; 
  inline double search_factor          () const { return _search_factor          ;} // =1.0; 
  inline double pT_min_leading_protojet() const { return _pT_min_leading_protojet;} // =0.; 
  inline double pT_min_second_protojet () const { return _pT_min_second_protojet ;} // =0.;
  inline int    merge_max              () const { return _merge_max              ;} // =10000; 
  inline double pT_min_nomerge         () const { return _pT_min_nomerge         ;} // =0.;


  /// access the split_ratio() also by the name overlap_threshold()
  inline double overlap_threshold() const {return split_ratio();}

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;
  /// the plugin mechanism's standard way of accessing the jet radius
  virtual double R() const {return cone_radius();}
  

private:

  double _cone_radius ;//= 0.5;
  double _min_jet_Et  ;//= 8.0;
  double _split_ratio ;//= 0.5; // overlap threshold
        
  //the parameters below have been found to be set to the values given below 
  //in the original implementation, shouldn't be altered
  double _far_ratio              ; // =0.5;
  double _Et_min_ratio           ; // =0.5;
  bool   _kill_duplicate         ; // =true;
  double _duplicate_dR           ; // =0.005; 
  double _duplicate_dPT          ; // =0.01; 
  double _search_factor          ; // =1.0; 
  double _pT_min_leading_protojet; // =0.; 
  double _pT_min_second_protojet ; // =0.;
  int    _merge_max              ; // =10000; 
  double _pT_min_nomerge         ; // =0.;

  // here are the variables for the default parameters of the D0 Run II Cone algorithm.
  // They are set in the .cc file 
  const static double _DEFAULT_split_ratio             ;// = 0.5  ; // overlap threshold
  const static double _DEFAULT_far_ratio               ;// = 0.5  ;
  const static double _DEFAULT_Et_min_ratio            ;// = 0.5  ;
  const static bool   _DEFAULT_kill_duplicate          ;// = true ;
  const static double _DEFAULT_duplicate_dR            ;// = 0.005; 
  const static double _DEFAULT_duplicate_dPT           ;// = 0.01 ; 
  const static double _DEFAULT_search_factor           ;// = 1.0  ; 
  const static double _DEFAULT_pT_min_leading_protojet ;// = 0.   ; 
  const static double _DEFAULT_pT_min_second_protojet  ;// = 0.   ;
  const static int    _DEFAULT_merge_max               ;// = 10000; 
  const static double _DEFAULT_pT_min_nomerge          ;// = 0.   ;

  static bool _first_time;

  /// print a banner for reference to the 3rd-party code
  void _print_banner(std::ostream *ostr) const;
};

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif // __D0RUNIICONEPLUGIN_HH__
