#ifndef __D0RUNIBASECONEPLUGIN_HH__
#define __D0RUNIBASECONEPLUGIN_HH__

//STARTHEADER
// $Id: D0RunIBaseConePlugin.hh 1778 2010-10-25 10:02:58Z soyez $
//
// Copyright (c) 2009-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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
/// @ingroup internal
/// \class D0RunIBaseConePlugin
///
/// D0RunIBaseConePlugin is base class for a plugin for FastJet (v3.0 or later) 
/// that provides an interface to the D0 version of Run-I cone algorithm
///
/// Note that this base class is purely virtual and thus needs to be
/// overloaded. In practice this means that you should use one of
/// D0RunIConePlugin or D0RunIpre96ConePlugin.
///
/// The D0 code has been obtained from Lars Sonnenschein's web-space
/// http://www-d0.fnal.gov/~sonne/D0RunIcone.tgz
///
/// The version of the D0 Run I code distributed here has been
/// modified by the FastJet authors, so as to provide access to the
/// contents of the jets (as is necessary for the plugin). This does
/// not modify the results of the clustering.
//
//----------------------------------------------------------------------
class D0RunIBaseConePlugin : public JetDefinition::Plugin {
public:
  /// A D0RunIConePlugin constructor which sets the "free" parameters of the
  /// algorithm:
  ///
  ///  \param CONErad is the cone radius
  ///
  ///  \param JETmne is a minimum ET requirement on every iteration
  ///    (jet dropped if Et < JETmne * Et_min_ratio ).
  ///    The value that has been used by D0 for JETmne: 8 GeV 
  ///    (and Et_min_ratio is 0.5)
  ///
  ///  \param SPlifr is the shared Et fraction splitting threshold, and
  ///    a value of 0.5 was usually used by D0
  ///
  /// The remaining parameters of the algorithm are not to be modified if the algorithm
  /// is to correspond to the one actually used by D0.
  D0RunIBaseConePlugin (double CONErad_in, 
			double JETmne_in, 
			double SPLifr_in = _DEFAULT_SPLifr) :
    _CONErad           (CONErad_in                 ),
    _JETmne            (JETmne_in                  ),
    _SPLifr            (SPLifr_in                  ),
    _TWOrad            (_DEFAULT_TWOrad            ),
    _D0_Angle          (_DEFAULT_D0_Angle          ),
    _Increase_Delta_R  (_DEFAULT_Increase_Delta_R  ),
    _Kill_Far_Clusters (_DEFAULT_Kill_Far_Clusters ),
    _Jet_Et_Min_On_Iter(_DEFAULT_Jet_Et_Min_On_Iter),
    _Far_Ratio         (_DEFAULT_Far_Ratio         ),
    _Eitem_Negdrop     (_DEFAULT_Eitem_Negdrop     ),
    _Et_Min_Ratio      (_DEFAULT_Et_Min_Ratio      ),
    _Thresh_Diff_Et    (_DEFAULT_Thresh_Diff_Et    ){}

  // some functions to return info about parameters
  inline double CONErad           () const { return _CONErad           ;} //= 0.7;
  inline double JETmne            () const { return _JETmne            ;} //= 8.0;
  inline double SPLifr            () const { return _SPLifr            ;} // =0.5;
  inline double TWOrad            () const { return _TWOrad            ;} //= 0.;
  inline bool   D0_Angle          () const { return _D0_Angle          ;} // =false;
  inline bool   Increase_Delta_R  () const { return _Increase_Delta_R  ;} // =true; 
  inline bool   Kill_Far_Clusters () const { return _Kill_Far_Clusters ;} // =true; 
  inline bool   Jet_Et_Min_On_Iter() const { return _Jet_Et_Min_On_Iter;} // =true; 
  inline double Far_Ratio         () const { return _Far_Ratio         ;} // =0.5; 
  inline double Eitem_Negdrop     () const { return _Eitem_Negdrop     ;} // =-1.0;
  inline double Et_Min_Ratio      () const { return _Et_Min_Ratio      ;} // =0.5; 
  inline double Thresh_Diff_Et    () const { return _Thresh_Diff_Et    ;} // =0.01;


  /// access the split_ratio() also by the name overlap_threshold()
  inline double overlap_threshold() const {return SPLifr();}

  // the things that are required by base class
  virtual std::string description () const = 0;

  // the part that really does the clustering
  virtual void run_clustering(ClusterSequence &) const = 0;

  /// the plugin mechanism's standard way of accessing the jet radius
  virtual double R() const {return CONErad();}

protected:
  template<typename HepEntityType>
  void run_clustering_worker(ClusterSequence &) const;

  //private:

  double _CONErad  ;//= 0.7 
  double _JETmne  ;//= 8.
  //the parameters below have been found to be set to the values given below 
  //in the original implementation, shouldn't be altered 
  double _SPLifr            ;  //=0.5
  double _TWOrad            ;  //=0.
  bool   _D0_Angle          ;  //=false
  bool   _Increase_Delta_R  ;  //=true
  bool   _Kill_Far_Clusters ;  //=true
  bool   _Jet_Et_Min_On_Iter;  //=true
  double _Far_Ratio         ;  //=0.5         
  double _Eitem_Negdrop     ;  //=-1.0
  double _Et_Min_Ratio      ;  //=0.5
  double _Thresh_Diff_Et    ;  //=0.01

  // here are the variables for the default parameters of the D0 Run I Cone algorithm.
  // They are set in the .cc file 
  const static double _DEFAULT_SPLifr                  ;  // = 0.5; //shared Et fraction threshold
  const static double _DEFAULT_TWOrad                  ;  // = 0.; //minimum Delta_R separation between cones
  const static bool   _DEFAULT_D0_Angle                ;  // = false;
  const static bool   _DEFAULT_Increase_Delta_R        ;  // = true;
  const static bool   _DEFAULT_Kill_Far_Clusters       ;  // = true;
  const static bool   _DEFAULT_Jet_Et_Min_On_Iter      ;  // = true;
  const static double _DEFAULT_Far_Ratio               ;  // = 0.5;
  const static double _DEFAULT_Eitem_Negdrop           ;  // = -1.0;
  const static double _DEFAULT_Et_Min_Ratio            ;  // = 0.5;
  const static double _DEFAULT_Thresh_Diff_Et          ;  // = 0.01;
};


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif // __D0RUNIBASECONEPLUGIN_HH__
