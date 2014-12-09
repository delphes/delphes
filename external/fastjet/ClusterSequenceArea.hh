//FJSTARTHEADER
// $Id: ClusterSequenceArea.hh 3484 2014-07-29 21:39:39Z soyez $
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

#ifndef __FASTJET_CLUSTERSEQUENCEAREA_HH__
#define __FASTJET_CLUSTERSEQUENCEAREA_HH__

#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceVoronoiArea.hh"
#include "fastjet/AreaDefinition.hh"

FASTJET_BEGIN_NAMESPACE

/// @ingroup area_classes
/// \class ClusterSequenceArea
/// General class for user to obtain ClusterSequence with additional
/// area information.
///
/// Based on the area_def, it automatically dispatches the work to the
/// appropriate actual ClusterSequenceAreaBase-derived-class to do the
/// real work.
class ClusterSequenceArea : public ClusterSequenceAreaBase {
public:
  /// main constructor
  template<class L> ClusterSequenceArea
         (const std::vector<L> & pseudojets, 
	  const JetDefinition & jet_def_in,
	  const AreaDefinition & area_def_in)  : _area_def(area_def_in) {
    initialize_and_run_cswa(pseudojets, jet_def_in);
  }

  /// constructor with a GhostedAreaSpec
  template<class L> ClusterSequenceArea
         (const std::vector<L> & pseudojets, 
	  const JetDefinition & jet_def_in,
	  const GhostedAreaSpec & ghost_spec)   : _area_def(ghost_spec){
    initialize_and_run_cswa(pseudojets, jet_def_in);
  }

  /// constructor with a VoronoiAreaSpec
  template<class L> ClusterSequenceArea
         (const std::vector<L> & pseudojets, 
	  const JetDefinition & jet_def_in,
	  const VoronoiAreaSpec & voronoi_spec)   : _area_def(voronoi_spec){
    initialize_and_run_cswa(pseudojets, jet_def_in);
  }

  /// return a reference to the area definition
  const AreaDefinition & area_def() const {return _area_def;}


  /// return the area associated with the given jet
  virtual double area       (const PseudoJet & jet) const {
    return _area_base->area(jet);}

  /// return the error (uncertainty) associated with the determination
  /// of the area of this jet
  virtual double area_error (const PseudoJet & jet) const {
    return _area_base->area_error(jet);}

  /// return the 4-vector area
  virtual PseudoJet area_4vector(const PseudoJet & jet) const {
    return _area_base->area_4vector(jet);}

  // /// return the total area, up to |y|<maxrap, that is free of jets
  // virtual double empty_area(double maxrap) const {
  //   return _area_base->empty_area(maxrap);}
  // 
  // /// return something similar to the number of pure ghost jets
  // /// in the given rapidity range in an active area case.
  // /// For the local implementation we return empty_area/(0.55 pi R^2),
  // /// based on measured properties of ghost jets with kt and cam. Note
  // /// that the number returned is a double.
  // virtual double n_empty_jets(double maxrap) const {
  //   return _area_base->n_empty_jets(maxrap);

  /// return the total area, corresponding to the given selector, that
  /// is free of jets
  ///
  /// The selector needs to have a finite area and be applicable jet by
  /// jet (see the BackgroundEstimator and Subtractor tools for more
  /// advanced usage)
  virtual double empty_area(const Selector & selector) const {
    return _area_base->empty_area(selector);}

  /// return something similar to the number of pure ghost jets
  /// in the given rap-phi range in an active area case.
  /// For the local implementation we return empty_area/(0.55 pi R^2),
  /// based on measured properties of ghost jets with kt and cam. Note
  /// that the number returned is a double.
  ///
  /// The selector needs to have a finite area and be applicable jet by
  /// jet (see the BackgroundEstimator and Subtractor tools for more
  /// advanced usage)
  virtual double n_empty_jets(const Selector & selector) const {
    return _area_base->n_empty_jets(selector);
  }

  /// true if a jet is made exclusively of ghosts
  virtual bool is_pure_ghost(const PseudoJet & jet) const {
    return _area_base->is_pure_ghost(jet);
  }

  /// true if this ClusterSequence has explicit ghosts
  virtual bool has_explicit_ghosts() const {
    return _area_base->has_explicit_ghosts();
  }
  

  /// overload version of what's in the ClusterSequenceAreaBase class, which 
  /// additionally checks compatibility between "selector" and region in which
  /// ghosts are thrown.
  ///
  /// The selector needs to have a finite area and be applicable jet by
  /// jet (see the BackgroundEstimator and Subtractor tools for more
  /// advanced usage)
  virtual void get_median_rho_and_sigma(const std::vector<PseudoJet> & all_jets,
					const Selector & selector, 
                                        bool use_area_4vector,
                                        double & median, double & sigma,
                                        double & mean_area,
					bool all_are_incl = false) const {
    _warn_if_range_unsuitable(selector);
    ClusterSequenceAreaBase::get_median_rho_and_sigma(
                                 all_jets, selector, use_area_4vector,
				 median, sigma, mean_area, all_are_incl);
  }

  /// overload version of what's in the ClusterSequenceAreaBase class,
  /// which actually just does the same thing as the base version (but
  /// since we've overridden the 5-argument version above, we have to
  /// override the 4-argument version too.
  virtual void get_median_rho_and_sigma(const Selector & selector, 
                                        bool use_area_4vector,
                                        double & median, double & sigma) const {
    ClusterSequenceAreaBase::get_median_rho_and_sigma(selector,use_area_4vector,
                                                      median,sigma);
  }

  /// overload version of what's in the ClusterSequenceAreaBase class,
  /// which actually just does the same thing as the base version (but
  /// since we've overridden the multi-argument version above, we have to
  /// override the 5-argument version too.
  virtual void get_median_rho_and_sigma(const Selector & selector, 
                                        bool use_area_4vector,
                                        double & median, double & sigma,
					double & mean_area) const {
    ClusterSequenceAreaBase::get_median_rho_and_sigma(selector,use_area_4vector,
                                                      median,sigma, mean_area);
  }


  /// overload version of what's in the ClusterSequenceAreaBase class, which 
  /// additionally checks compatibility between "range" and region in which
  /// ghosts are thrown.
  virtual void parabolic_pt_per_unit_area(double & a, double & b, 
                                          const Selector & selector, 
                                          double exclude_above=-1.0, 
                                          bool use_area_4vector=false) const {
    _warn_if_range_unsuitable(selector);
    ClusterSequenceAreaBase::parabolic_pt_per_unit_area(
                                a,b,selector, exclude_above, use_area_4vector);
  }


private:
  
  /// print a warning if the range is unsuitable for the current
  /// calculation of the area (e.g. because ghosts do not extend
  /// far enough).
  void _warn_if_range_unsuitable(const Selector & selector) const;

  template<class L> void initialize_and_run_cswa (
                                 const std::vector<L> & pseudojets, 
                                 const JetDefinition & jet_def);

  std::auto_ptr<ClusterSequenceAreaBase> _area_base;
  AreaDefinition _area_def;
  static LimitedWarning _range_warnings;
  static LimitedWarning _explicit_ghosts_repeats_warnings;

};

//----------------------------------------------------------------------
template<class L> void ClusterSequenceArea::initialize_and_run_cswa(
           const std::vector<L> & pseudojets, 
           const JetDefinition  & jet_def_in)
 {
  
  ClusterSequenceAreaBase * _area_base_ptr;
  switch(_area_def.area_type()) {
  case active_area:
    _area_base_ptr = new ClusterSequenceActiveArea(pseudojets, 
                                                   jet_def_in, 
                                                   _area_def.ghost_spec());
    break;
  case active_area_explicit_ghosts:
    if (_area_def.ghost_spec().repeat() != 1) 
      _explicit_ghosts_repeats_warnings.warn("Requested active area with explicit ghosts with repeat != 1; only 1 set of ghosts will be used");
    _area_base_ptr = new ClusterSequenceActiveAreaExplicitGhosts(pseudojets, 
                                                   jet_def_in, 
                                                   _area_def.ghost_spec());
    break;
  case voronoi_area:
    _area_base_ptr = new ClusterSequenceVoronoiArea(pseudojets, 
                                                   jet_def_in, 
                                                   _area_def.voronoi_spec());
    break;
  case one_ghost_passive_area:
    _area_base_ptr = new ClusterSequence1GhostPassiveArea(pseudojets, 
						    jet_def_in, 
						    _area_def.ghost_spec());
    break;
  case passive_area:
    _area_base_ptr = new ClusterSequencePassiveArea(pseudojets, 
						    jet_def_in, 
						    _area_def.ghost_spec());
    break;
  default:
    std::ostringstream err;
    err << "Error: unrecognized area_type in ClusterSequenceArea:" 
	<< _area_def.area_type();
    throw Error(err.str());
    //exit(-1);
  }
  // now copy across the information from the area base class
  _area_base = std::auto_ptr<ClusterSequenceAreaBase>(_area_base_ptr);
  transfer_from_sequence(*_area_base);
}

FASTJET_END_NAMESPACE

#endif // __FASTJET_CLUSTERSEQUENCEAREA_HH__


