#ifndef __FASTJET_BACKGROUND_ESTIMATOR_HH__
#define __FASTJET_BACKGROUND_ESTIMATOR_HH__


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

#include "fastjet/config.h"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/BackgroundEstimatorBase.hh"
#include <iostream>

#ifdef FASTJET_HAVE_THREAD_SAFETY
#include <atomic>
#endif

FASTJET_BEGIN_NAMESPACE     // defined in fastjet/internal/base.hh


/// @ingroup tools_background
/// \class JetMedianBackgroundEstimator
///
/// Class to estimate the pt density of the background per unit area,
/// using the median of the distribution of pt/area from jets that
/// pass some selection criterion.
///
/// Events are passed either in the form of the event particles (in
/// which they're clustered by the class), a ClusterSequenceArea (in
/// which case the jets used are those returned by "inclusive_jets()")
/// or directly as a set of jets.
///
/// The selection criterion is typically a geometrical one (e.g. all
/// jets with |y|<2) sometimes supplemented with some kinematical
/// restriction (e.g. exclusion of the two hardest jets). It is passed
/// to the class through a Selector.
///
/// Beware: 
///   by default, to correctly handle partially empty events, the
///   class attempts to calculate an "empty area", based
///   (schematically) on
///
///          range.total_area() - sum_{jets_in_range} jets.area()
///  
///   For ranges with small areas, this can be inaccurate (particularly 
///   relevant in dense events where empty_area should be zero and ends
///   up not being zero).
///
///   This calculation of empty area can be avoided if a
///   ClusterSequenceArea class with explicit ghosts
///   (ActiveAreaExplicitGhosts) is used.  This is _recommended_
///   unless speed requirements cause you to use Voronoi areas. For
///   speedy background estimation you could also consider using
///   GridMedianBackgroundEstimator.
///
///
class JetMedianBackgroundEstimator : public BackgroundEstimatorBase {
public:
  /// @name constructors and destructors
  //\{
  //----------------------------------------------------------------
  /// Constructor that sets the rho range as well as the jet
  /// definition and area definition to be used to cluster the
  /// particles. Prior to the estimation of rho, one has to provide
  /// the particles to cluster using set_particles(...)
  ///
  /// \param rho_range  the Selector specifying which jets will be considered
  /// \param jet_def    the jet definition to use for the clustering
  /// \param area_def   the area definition to use for the clustering
  JetMedianBackgroundEstimator(const Selector &rho_range,
                               const JetDefinition &jet_def,
                               const AreaDefinition &area_def);

  /// ctor from a ClusterSequenceAreaBase with area
  ///
  /// \param rho_range   the Selector specifying which jets will be considered
  /// \param csa         the ClusterSequenceArea to use
  ///
  /// Pre-conditions: 
  ///  - one should be able to estimate the "empty area" (i.e. the area
  ///    not occupied by jets). This is feasible if at least one of the following
  ///    conditions is satisfied:
  ///     ( i) the ClusterSequence has explicit ghosts
  ///     (ii) the range has a computable area.
  ///  - the jet algorithm must be suited for median computation
  ///    (otherwise a warning will be issues)
  ///
  /// Note that selectors with e.g. hardest-jets exclusion do not have
  /// a well-defined area. For this reasons, it is STRONGLY advised to
  /// use an area with explicit ghosts.
  JetMedianBackgroundEstimator(const Selector &rho_range,
                               const ClusterSequenceAreaBase &csa);


  /// Default constructor that optionally sets the rho range. The
  /// configuration must be done later calling
  /// set_cluster_sequence(...) or set_jets(...).
  ///
  /// \param rho_range   the Selector specifying which jets will be considered
  ///
  JetMedianBackgroundEstimator(const Selector &rho_range = SelectorIdentity())
    : _rho_range(rho_range), _jet_def(JetDefinition()),
      _enable_rho_m(true){ reset(); }

  
  /// default dtor
  ~JetMedianBackgroundEstimator(){}

  //\}


  /// @name setting a new event
  //\{
  //----------------------------------------------------------------

  /// tell the background estimator that it has a new event, composed
  /// of the specified particles.
  virtual void set_particles(const std::vector<PseudoJet> & particles) FASTJET_OVERRIDE;

  // tell the background estimator that it has a new event, composed
  // of the specified particles and use the supplied seed for the
  // generation of ghosts. If the seed is empty, it is ignored.
  virtual void set_particles_with_seed(const std::vector<PseudoJet> & particles, const std::vector<int> & seed) FASTJET_OVERRIDE;

  /// (re)set the cluster sequence (with area support) to be used by
  /// future calls to rho() etc. 
  ///
  /// \param csa  the cluster sequence area
  ///
  /// Pre-conditions: 
  ///  - one should be able to estimate the "empty area" (i.e. the area
  ///    not occupied by jets). This is feasible if at least one of the following
  ///    conditions is satisfied:
  ///     ( i) the ClusterSequence has explicit ghosts
  ///     (ii) the range selected has a computable area.
  ///  - the jet algorithm must be suited for median computation
  ///    (otherwise a warning will be issues)
  ///
  /// Note that selectors with e.g. hardest-jets exclusion do not have
  /// a well-defined area. For this reasons, it is STRONGLY advised to
  /// use an area with explicit ghosts.
  void set_cluster_sequence(const ClusterSequenceAreaBase & csa);

  /// (re)set the jets (which must have area support) to be used by future
  /// calls to rho() etc.; for the conditions that must be satisfied
  /// by the jets, see the Constructor that takes jets.
  void set_jets(const std::vector<PseudoJet> &jets);

  /// (re)set the selector to be used for future calls to rho() etc.
  void set_selector(const Selector & rho_range_selector) {
    _rho_range = rho_range_selector;
    _set_cache_unavailable();
  }

  /// determine whether the automatic calculation of rho_m and sigma_m
  /// is enabled (by default true)
  void set_compute_rho_m(bool enable){
    _enable_rho_m = enable;
    _set_cache_unavailable();
  }

  //\}


  /// return a pointer to a copy of this BGE; the user is responsible
  /// for eventually deleting the resulting object.
  BackgroundEstimatorBase * copy() const FASTJET_OVERRIDE {
    return new JetMedianBackgroundEstimator(*this);
  };


  /// @name  retrieving fundamental information
  //\{
  //----------------------------------------------------------------
  /// get the full set of background properties
  ///
  /// For background estimators using a local ranges, this throws an
  ///   error (use estimate(jet) instead)
  /// In the presence of a rescaling, the rescaling is not included
  BackgroundEstimate estimate() const FASTJET_OVERRIDE;
  
  /// get the full set of background properties for a given reference jet
  BackgroundEstimate estimate(const PseudoJet &jet) const FASTJET_OVERRIDE;

  /// get rho, the median background density per unit area
  double rho() const FASTJET_OVERRIDE;

  /// get sigma, the background fluctuations per unit area
  double sigma() const FASTJET_OVERRIDE;

  /// get rho, the median background density per unit area, locally at
  /// the position of a given jet.
  ///
  /// If the Selector associated with the range takes a reference jet
  /// (i.e. is relocatable), then for subsequent operations the
  /// Selector has that jet set as its reference.
  double rho(const PseudoJet & jet) FASTJET_OVERRIDE;

  /// get sigma, the background fluctuations per unit area,
  /// locally at the position of a given jet.
  ///
  /// If the Selector associated with the range takes a reference jet
  /// (i.e. is relocatable), then for subsequent operations the
  /// Selector has that jet set as its reference.
  double sigma(const PseudoJet &jet) FASTJET_OVERRIDE;

  /// returns true if this background estimator has support for
  /// determination of sigma
  virtual bool has_sigma() const FASTJET_OVERRIDE {return true;}

  //----------------------------------------------------------------
  // now do the same thing for rho_m and sigma_m

  /// returns rho_m, the purely longitudinal, particle-mass-induced
  /// component of the background density per unit area
  virtual double rho_m() const FASTJET_OVERRIDE;

  /// returns sigma_m, a measure of the fluctuations in the purely
  /// longitudinal, particle-mass-induced component of the background
  /// density per unit area; must be multipled by sqrt(area) to get
  /// fluctuations for a region of a given area.
  virtual double sigma_m() const FASTJET_OVERRIDE;

  /// Returns rho_m locally at the jet position. As for rho(jet), it is non-const.
  virtual double rho_m(const PseudoJet & /*jet*/) FASTJET_OVERRIDE;

  /// Returns sigma_m locally at the jet position. As for rho(jet), it is non-const.
  virtual double sigma_m(const PseudoJet & /*jet*/) FASTJET_OVERRIDE;

  /// Returns true if this background estimator has support for
  /// determination of rho_m.
  ///
  /// In te presence of a density class, support for rho_m is
  /// automatically disabled
  ///
  /// Note that support for sigma_m is automatic is one has sigma and
  /// rho_m support.
  virtual bool has_rho_m() const FASTJET_OVERRIDE {return _enable_rho_m && (_jet_density_class == 0);}
  //\}
  
  /// @name  retrieving additional useful information
  //\{
  //----------------------------------------------------------------
  /// Returns the mean area of the jets used to actually compute the
  /// background properties in the last call of rho() or sigma()
  /// If the configuration has changed in the meantime, throw an error.
  double mean_area() const;
  
  /// returns the number of jets used to actually compute the
  /// background properties in the last call of rho() or sigma()
  /// If the configuration has changed in the meantime, throw an error.
  unsigned int n_jets_used() const;

  /// returns the jets used to actually compute the background
  /// properties
  std::vector<PseudoJet> jets_used() const;

  /// Returns the estimate of the area (within the range defined by
  /// the selector) that is not occupied by jets. The value is that
  /// for the last call of rho() or sigma()
  /// If the configuration has changed in the meantime, throw an error.
  ///
  /// The answer is defined to be zero if the area calculation
  /// involved explicit ghosts; if the area calculation was an active
  /// area, then use is made of the active area's internal list of
  /// pure ghost jets (taking those that pass the selector); otherwise
  /// it is based on the difference between the selector's total area
  /// and the area of the jets that pass the selector.
  ///
  /// The result here is just the cached result of the corresponding
  /// call to the ClusterSequenceAreaBase function.
  double empty_area() const;

  /// Returns the number of empty jets used when computing the
  /// background properties. The value is that for the last call of
  /// rho() or sigma().
  /// If the configuration has changed in the meantime, throw an error.
  ///
  /// If the area has explicit ghosts the result is zero; for active
  /// areas it is the number of internal pure ghost jets that pass the
  /// selector; otherwise it is deduced from the empty area, divided by 
  /// \f$ 0.55 \pi R^2 \f$ (the average pure-ghost-jet area).
  ///
  /// The result here is just the cached result of the corresponding
  /// call to the ClusterSequenceAreaBase function.
  double n_empty_jets() const;

  //}


  /// @name configuring behaviour
  //\{
  //----------------------------------------------------------------

  /// Resets the class to its default state, including the choice to
  /// use 4-vector areas.
  ///
  void reset();

  /// By default when calculating pt/Area for a jet, it is the
  /// transverse component of the 4-vector area that is used in the ratiof \f$p_t/A\f$. 
  /// Calling this function with a "false" argument causes the scalar area to
  /// be used instead. 
  ///
  /// While the difference between the two choices is usually small,
  /// for high-precision work it is usually the 4-vector area that is
  /// to be preferred.
  ///
  ///  \param use_it             whether one uses the 4-vector area or not (true by default)
  void set_use_area_4vector(bool use_it = true){
    _use_area_4vector = use_it;
    _set_cache_unavailable();
  }  

  /// check if the estimator uses the 4-vector area or the scalar area
  bool use_area_4vector() const{ return _use_area_4vector;}

  /// The FastJet v2.X sigma calculation had a small spurious offset
  /// in the limit of a small number of jets. This is fixed by default
  /// in versions 3 upwards. The old behaviour can be obtained with a
  /// call to this function.
  void set_provide_fj2_sigma(bool provide_fj2_sigma = true) {
    _provide_fj2_sigma = provide_fj2_sigma;
    _set_cache_unavailable();
  }

  /// Set a pointer to a class that calculates the quantity whose
  /// median will be calculated; if the pointer is null then pt/area
  /// is used (as occurs also if this function is not called).
  ///
  /// Note that this is still <i>preliminary</i> in FastJet 3.0 and
  /// that backward compatibility is not guaranteed in future releases
  /// of FastJet
  void set_jet_density_class(const FunctionOfPseudoJet<double> * jet_density_class);

  /// return the pointer to the jet density class
  const FunctionOfPseudoJet<double> *  jet_density_class() const{
    return _jet_density_class;
  }

  /// Set a pointer to a class that calculates the rescaling factor as
  /// a function of the jet (position). Note that the rescaling factor
  /// is used both in the determination of the "global" rho (the pt/A
  /// of each jet is divided by this factor) and when asking for a
  /// local rho (the result is multiplied by this factor).
  ///
  /// The BackgroundRescalingYPolynomial class can be used to get a
  /// rescaling that depends just on rapidity.
  virtual void set_rescaling_class(const FunctionOfPseudoJet<double> * rescaling_class_in) FASTJET_OVERRIDE {
    BackgroundEstimatorBase::set_rescaling_class(rescaling_class_in);
    _set_cache_unavailable();
  }

  //\}

  /// @name description
  //\{
  //----------------------------------------------------------------

  /// returns a textual description of the background estimator
  std::string description() const FASTJET_OVERRIDE;

  //\}

  /// an internal class to hold results of the calculation
  /// that are to be assigned to the "extras" part of a BackgroundEstimate
  class Extras : public BackgroundEstimate::Extras {
  public:
    Extras()
      :  _reference_jet(PseudoJet()), _n_jets_used(0),
         _n_empty_jets(0.0), _empty_area(0.0) {}

    /// returns the current reference jet
    PseudoJet reference_jet() const {return _reference_jet;}
    
    /// returns the number of jets used to estimate the background
    unsigned int n_jets_used() const {return _n_jets_used;}
    
    /// returns the number of empty (pure-ghost) jets
    double n_empty_jets() const {return _n_empty_jets;}
      
    /// returns the empty (pure-ghost/unclustered) area!
    double empty_area() const {return _empty_area;}

    void set_reference_jet(const PseudoJet &reference_jet_in){
      _reference_jet = reference_jet_in;
    }
    void set_n_jets_used(int n_jets_used_in){ _n_jets_used=n_jets_used_in;}
    void set_n_empty_jets(double n_empty_jets_in){ _n_empty_jets=n_empty_jets_in;}
    void set_empty_area(double empty_area_in){ _empty_area=empty_area_in;}
    
    
  protected:
    PseudoJet _reference_jet;  ///< current reference jet
    unsigned int _n_jets_used; ///< number of jets used to estimate the background
    double _n_empty_jets;      ///< number of empty (pure-ghost) jets
    double _empty_area;        ///< the empty (pure-ghost/unclustered) area!
  };

  
private:

  /// compute the background properties for a given jet (excluding
  /// rescaling factors) and return a corresponding BackgroundEstimate
  ///
  /// this leaves the cache (and the status flags) unchanged
  BackgroundEstimate _compute(const PseudoJet &jet) const;

  //------
  // the next calls are meant for the case where the cache can be
  // filled once and for all, i.e. cases where the selector does NOT
  // take a reference
  
  /// fill the cache with the given estimate
  void _cache_no_overwrite(const BackgroundEstimate &estimate) const;
   
  /// fill the cache with a computed estimate
  void _compute_and_cache_no_overwrite() const;

  //------
  // the next calls are meant for the case where the selector does
  // take a reference and the cache needs to be refilled whenever one
  // calls this background estimate with a different reference jet

  /// fill the cache with the given estimate
  void _cache(const BackgroundEstimate &estimate) const;
   
  /// update the cache if need be and return the background
  /// estimate. This is meant to be called for cases with a local
  /// range (selector that takesa reference)
  BackgroundEstimate _compute_and_cache_if_needed(const PseudoJet &jet) const;
  //-----
  
  /// check that the underlying structure is still alive
  /// throw an error otherwise
  void _check_csa_alive() const;

  /// check that the algorithm used for the clustering is adapted for
  /// background estimation (i.e. either kt or C/A)
  /// Issue a warning otherwise
  void _check_jet_alg_good_for_median() const;

  // the basic parameters of this class (passed through the variou ctors)
  Selector _rho_range;                   ///< range to compute the background in
  JetDefinition _jet_def;                ///< the jet def to use for teh clustering
  AreaDefinition _area_def;              ///< the area def to use for teh clustering
  std::vector<PseudoJet> _included_jets; ///< jets to be used
  
  // the tunable parameters of the class
  bool _use_area_4vector;
  bool _provide_fj2_sigma;
  const FunctionOfPseudoJet<double> * _jet_density_class;
  //SharedPtr<BackgroundRescalingBase> _rescaling_class_sharedptr;
  bool _enable_rho_m;
  
  // internal variables
  SharedPtr<PseudoJetStructureBase> _csi; ///< allows to check if _csa is still valid

  /// handle warning messages
  static LimitedWarning _warnings;
  static LimitedWarning _warnings_zero_area;
  static LimitedWarning _warnings_preliminary;
};




//----------------------------------------------------------------------
/// @ingroup tools_background
/// \class BackgroundJetPtDensity
/// Class that implements pt/area_4vector.perp() for background estimation
/// <i>(this is a preliminary class)</i>.
class BackgroundJetPtDensity : public FunctionOfPseudoJet<double> {
public:
  virtual double result(const PseudoJet & jet) const {
    return jet.perp() / jet.area_4vector().perp();
  }
  virtual std::string description() const {return "BackgroundJetPtDensity";}
};


//----------------------------------------------------------------------
/// @ingroup tools_background
/// \class BackgroundJetScalarPtDensity
/// Class that implements (scalar pt sum of jet)/(scalar area of jet)
/// for background estimation <i>(this is a preliminary class)</i>.
///
/// Optionally it can return a quantity based on the sum of pt^n,
/// e.g. for use in subtracting fragementation function moments.
class BackgroundJetScalarPtDensity : public FunctionOfPseudoJet<double> {
public:
  /// Default constructor provides background estimation with scalar pt sum
  BackgroundJetScalarPtDensity() : _pt_power(1) {}

  /// Constructor to provide background estimation based on 
  /// \f$ sum_{i\in jet} p_{ti}^{n} \f$
  BackgroundJetScalarPtDensity(double n) : _pt_power(n) {}

  virtual double result(const PseudoJet & jet) const;

  virtual std::string description() const;

private:
  double _pt_power;
};

//----------------------------------------------------------------------
/// @ingroup tools_background
/// \class BackgroundJetPtMDensity
/// Class that implements
/// \f$  \frac{1}{A} \sum_{i \in jet} (\sqrt{p_{ti}^2+m^2} - p_{ti}) \f$
/// for background estimation <i>(this is a preliminary class)</i>.
/// 
///
/// This is useful for correcting jet masses in cases where the event
/// involves massive particles.
class BackgroundJetPtMDensity : public FunctionOfPseudoJet<double> {
public:
  virtual double result(const PseudoJet & jet) const {
    std::vector<PseudoJet> constituents = jet.constituents();
    double scalar_ptm = 0;
    for (unsigned i = 0; i < constituents.size(); i++) {
      scalar_ptm += constituents[i].mperp() - constituents[i].perp();
    }
    return scalar_ptm / jet.area();
  }

  virtual std::string description() const {return "BackgroundPtMDensity";}
};


FASTJET_END_NAMESPACE

#endif  // __BACKGROUND_ESTIMATOR_HH__

