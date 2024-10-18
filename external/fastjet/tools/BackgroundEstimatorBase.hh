#ifndef __FASTJET_BACKGROUND_ESTIMATOR_BASE_HH__
#define __FASTJET_BACKGROUND_ESTIMATOR_BASE_HH__

//FJSTARTHEADER
// $Id$
//
// Copyright (c) 2005-2024, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/Error.hh"
#include <iostream>

FASTJET_BEGIN_NAMESPACE     // defined in fastjet/internal/base.hh


/// @ingroup tools_background
/// @name helpers to handle the result of the background estimation
//\{
///
/// /// a class that holds the result of the calculation
///
/// By default it provides access to the main background properties:
/// rho, rho_m, sigma and sigma_m. If background estimators derived
/// from the base class want to store more information, this can be
/// done using the "Extra" information.
class BackgroundEstimate{
public:
  /// ctor wo initialisation
  BackgroundEstimate()
    : _rho(0.0), _sigma(0.0), _rho_m(0.0), _sigma_m(0.0), 
      _has_sigma(false), _has_rho_m(false),
      _mean_area(0.0){}

  
  /// @name for accessing information about the background
  ///@{

  /// background density per unit area
  double rho() const {return _rho;}

  /// background fluctuations per unit square-root area
  /// must be multipled by sqrt(area) to get fluctuations for a region
  /// of a given area.
  double sigma() const {return _sigma;}

  /// true if this background estimate has a determination of sigma
  bool has_sigma() {return true;}

  /// purely longitudinal (particle-mass-induced)
  /// component of the background density per unit area
  double rho_m() const {return _rho_m;}

  /// fluctuations in the purely longitudinal (particle-mass-induced)
  /// component of the background density per unit square-root area
  double sigma_m() const {return _sigma_m;}

  /// true if this background estimate has a determination of rho_m.
  /// Support for sigma_m is automatic if one has sigma and rho_m support.
  bool has_rho_m() const {return _has_rho_m;}

  /// mean area of the patches used to compute the background properties
  double mean_area() const {return _mean_area;}

  /// base class for extra information
  class Extras {
  public:
    // dummy ctor
    Extras(){};

    // dummy virtual dtor
    // makes it polymorphic to allow for dynamic_cast
    virtual ~Extras(){}; 
  };

  /// returns true if the background estimate has extra info
  bool has_extras() const{
    return _extras.get();
  }
  
  /// returns true if the background estimate has extra info
  /// compatible with the provided template type
  template<typename T>
  bool has_extras() const{
    return _extras.get() && dynamic_cast<const T *>(_extras.get());
  }

  /// returns a reference to the extra information associated with a
  /// given BackgroundEstimator. It assumes that the extra
  /// information is reachable with class name
  /// BackgroundEstimator::Extras
  template<typename BackgroundEstimator>
  const typename BackgroundEstimator::Extras & extras() const{
    return dynamic_cast<const typename BackgroundEstimator::Extras &>(* _extras.get());
  }

  ///@}


  /// @name for setting information about the background (internal FJ use)
  ///@{

  /// reset to default
  void reset(){
    _rho = _sigma = _rho_m = _sigma_m = _mean_area = 0.0;
    _has_sigma = _has_rho_m = false;
    _extras.reset();
  }
  void set_rho(double rho_in) {_rho = rho_in;}
  void set_sigma(double sigma_in) {_sigma = sigma_in;}
  void set_has_sigma(bool has_sigma_in) {_has_sigma = has_sigma_in;}
  void set_rho_m(double rho_m_in) {_rho_m = rho_m_in;}
  void set_sigma_m(double sigma_m_in) {_sigma_m = sigma_m_in;}
  void set_has_rho_m(bool has_rho_m_in) {_has_rho_m = has_rho_m_in;}
  void set_mean_area(double mean_area_in) {_mean_area = mean_area_in;}

  /// apply a rescaling factor (to rho, rho_m, sigma, sigma_m)
  void apply_rescaling_factor(double rescaling_factor){
    _rho     *= rescaling_factor;
    _sigma   *= rescaling_factor;
    _rho_m   *= rescaling_factor;
    _sigma_m *= rescaling_factor;
  }

  /// sets the extra info based on the provided pointer
  ///
  /// When calling this method, the BackgroundEstimate class takes
  /// ownership of the pointer (and is responsible for deleting it)
  void set_extras(Extras *extras_in) {
    _extras.reset(extras_in);
  }
  ///@}


protected:
  double _rho;       ///< background estimated density per unit area
  double _sigma;     ///< background estimated fluctuations
  double _rho_m;     ///< "mass" background estimated density per unit area
  double _sigma_m;   ///< "mass" background estimated fluctuations
  bool _has_sigma;   ///< true if this estimate has a determination of sigma
  bool _has_rho_m;   ///< true if this estimate has a determination of rho_m
  double _mean_area; ///< mean area of the patches used to compute the bkg properties
  

  SharedPtr<Extras> _extras;

};


/// @ingroup tools_background
/// \class BackgroundEstimatorBase
///
/// Abstract base class that provides the basic interface for classes
/// that estimate levels of background radiation in hadron and
/// heavy-ion collider events.
///
class BackgroundEstimatorBase {
public:
  /// @name  constructors and destructors
  //\{
  //----------------------------------------------------------------
  BackgroundEstimatorBase() : _rescaling_class(0) {
    _set_cache_unavailable();                                               
  }

#ifdef FASTJET_HAVE_THREAD_SAFETY
  /// because of the internal atomic variale, we need to explicitly
  /// implement a copy ctor
  BackgroundEstimatorBase(const BackgroundEstimatorBase &other_bge);
#endif

  /// a default virtual destructor that does nothing
  virtual ~BackgroundEstimatorBase() {}
  //\}

  /// @name setting a new event
  //\{
  //----------------------------------------------------------------

  /// tell the background estimator that it has a new event, composed
  /// of the specified particles.
  virtual void set_particles(const std::vector<PseudoJet> & particles) = 0;

  /// an alternative call that takes a random number generator seed
  /// (typically a vector of length 2) to ensure reproducibility of
  /// background estimators that rely on random numbers (specifically
  /// JetMedianBackgroundEstimator with ghosted areas)
  virtual void set_particles_with_seed(const std::vector<PseudoJet> & particles, const std::vector<int> & /*seed*/) {
    set_particles(particles);
  }

  //\}

  /// return a pointer to a copy of this BGE; the user is responsible
  /// for eventually deleting the resulting object.
  virtual BackgroundEstimatorBase * copy() const = 0;

  /// @name  retrieving fundamental information
  //\{
  //----------------------------------------------------------------
  /// get the full set of background properties
  virtual BackgroundEstimate estimate() const = 0;
  
  /// get the full set of background properties for a given reference jet
  virtual BackgroundEstimate estimate(const PseudoJet &jet) const = 0;

  /// get rho, the background density per unit area
  virtual double rho() const = 0;

  /// get sigma, the background fluctuations per unit square-root area;
  /// must be multipled by sqrt(area) to get fluctuations for a region
  /// of a given area.
  virtual double sigma() const { 
    throw Error("sigma() not supported for this Background Estimator");
  }

  /// get rho, the background density per unit area, locally at the
  /// position of a given jet. Note that this is not const, because a
  /// user may then wish to query other aspects of the background that
  /// could depend on the position of the jet last used for a rho(jet)
  /// determination.
  virtual double rho(const PseudoJet & jet) = 0;

  /// get sigma, the background fluctuations per unit area, locally at
  /// the position of a given jet. As for rho(jet), it is non-const.
  virtual double sigma(const PseudoJet & /*jet*/) { 
    throw Error("sigma(jet) not supported for this Background Estimator");
  }

  /// returns true if this background estimator has support for
  /// determination of sigma
  virtual bool has_sigma() const {return false;}

  //----------------------------------------------------------------
  // now do the same thing for rho_m and sigma_m

  /// returns rho_m, the purely longitudinal, particle-mass-induced
  /// component of the background density per unit area
  virtual double rho_m() const{
    throw Error("rho_m() not supported for this Background Estimator");
  }

  /// returns sigma_m, a measure of the fluctuations in the purely
  /// longitudinal, particle-mass-induced component of the background
  /// density per unit square-root area; must be multipled by sqrt(area) to get
  /// fluctuations for a region of a given area.
  virtual double sigma_m() const { 
    throw Error("sigma_m() not supported for this Background Estimator");
  }

  /// Returns rho_m locally at the jet position. As for rho(jet), it is non-const.
  virtual double rho_m(const PseudoJet & /*jet*/){
    throw Error("rho_m(jet) not supported for this Background Estimator");
  }

  /// Returns sigma_m locally at the jet position. As for rho(jet), it is non-const.
  virtual double sigma_m(const PseudoJet & /*jet*/) { 
    throw Error("sigma_m(jet) not supported for this Background Estimator");
  }

  /// Returns true if this background estimator has support for
  /// determination of rho_m.
  ///
  /// Note that support for sigma_m is automatic is one has sigma and
  /// rho_m support.
  virtual bool has_rho_m() const {return false;}
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
  /// There is currently no support for different rescaling classes 
  /// for rho and rho_m determinations.
  virtual void set_rescaling_class(const FunctionOfPseudoJet<double> * rescaling_class_in) { _rescaling_class = rescaling_class_in; }

  /// return the pointer to the jet density class
  const FunctionOfPseudoJet<double> *  rescaling_class() const{
    return _rescaling_class;
  }

  //\}

  /// @name description
  //\{
  //----------------------------------------------------------------

  /// returns a textual description of the background estimator
  virtual std::string description() const = 0;

  //\}

  
protected:
  /// @name helpers for derived classes
  ///
  /// Note that these helpers are related to median-based estimation
  /// of the background, so there is no guarantee that they will
  /// remain in this base class in the long term
  //\{
  //----------------------------------------------------------------

  /// given a quantity in a vector (e.g. pt_over_area) and knowledge
  /// about the number of empty jets, calculate the median and
  /// stand_dev_if_gaussian (roughly from the 16th percentile)
  ///
  /// If do_fj2_calculation is set to true then this performs FastJet
  /// 2.X estimation of the standard deviation, which has a spurious
  /// offset in the limit of a small number of jets.
  void _median_and_stddev(const std::vector<double> & quantity_vector, 
                          double n_empty_jets, 
                          double & median, 
                          double & stand_dev_if_gaussian,
                          bool do_fj2_calculation = false
                          ) const;

  /// computes a percentile of a given _sorted_ vector
  ///  \param sorted_quantity_vector   the vector contains the data sample
  ///  \param percentile               the percentile (defined between 0 and 1) to compute
  ///  \param nempty                   an additional number of 0's
  ///                                  (considered at the beginning of 
  ///                                  the quantity vector)
  ///  \param do_fj2_calculation       carry out the calculation as it
  ///                                  was done in fj2 (suffers from "edge effects")
  double _percentile(const std::vector<double> & sorted_quantity_vector, 
                     const double percentile, 
                     const double nempty=0.0,
                     const bool do_fj2_calculation = false) const;

  //\}

  const FunctionOfPseudoJet<double> * _rescaling_class;
  static LimitedWarning _warnings_empty_area;


  // cached actual results of the computation
  mutable bool _cache_available;
  mutable BackgroundEstimate _cached_estimate;    ///< all the info about what is computed
  //\}

  /// @name helpers for thread safety
  ///
  /// Note that these helpers are related to median-based estimation
  /// of the background, so there is no guarantee that they will
  /// remain in this base class in the long term
  //\{
#ifdef FASTJET_HAVE_THREAD_SAFETY
  // true when one is currently writing to cache (i.e. when the spin lock is set)
  mutable std::atomic<bool> _writing_to_cache;

  void _set_cache_unavailable(){
    _cache_available = false;
    _writing_to_cache = false;
  }

  // // allows us to lock things down before caching basic (patches) info  
  // std::mutex _jets_caching_mutex;
#else
  void _set_cache_unavailable(){
    _cache_available = false;
  }
#endif

  void _lock_if_needed()   const;
  void _unlock_if_needed() const;

  //\}

};



//----------------------------------------------------------------------
/// @ingroup tools_background
/// A background rescaling that is a simple polynomial in y
class BackgroundRescalingYPolynomial : public FunctionOfPseudoJet<double> {
public:
  /// construct a background rescaling polynomial of the form
  /// a0 + a1*y + a2*y^2 + a3*y^3 + a4*y^4
  ///
  /// The following values give a reasonable reproduction of the
  /// Pythia8 tune 4C background shape for pp collisions at
  /// sqrt(s)=7TeV:
  ///
  /// - a0 =  1.157
  /// - a1 =  0
  /// - a2 = -0.0266
  /// - a3 =  0
  /// - a4 =  0.000048
  ///
  BackgroundRescalingYPolynomial(double a0=1, 
                                 double a1=0, 
                                 double a2=0, 
                                 double a3=0, 
                                 double a4=0) : _a0(a0), _a1(a1), _a2(a2), _a3(a3), _a4(a4) {}

  /// return the rescaling factor associated with this jet
  virtual double result(const PseudoJet & jet) const;
private:
  double _a0, _a1, _a2, _a3, _a4;
};





FASTJET_END_NAMESPACE

#endif  // __BACKGROUND_ESTIMATOR_BASE_HH__

// //--------------------------------------------------
// // deprecated
// class JetMedianBGE{
//   BackgroundEstimateDefinition();
// 
//   ....;
// }
// //--------------------------------------------------
// 
// class BackgroundEstimateDefinition{ 
//   //const EventBackground get_event_background(particles, <seed>) const;
// 
// 
//   //--------------------------------------------------
//   // DEPRECATED
//   void set_particles() {
// 
//     _worker = ...;
//   double rho(const PseudoJet &/*jet*/) const{ _worker.rho();}
// 
//   //--------------------------------------------------
//   
//   EventBackground(Worker?) _cached_event_background;
// };
// 
// class EventBackground{
//   EventBackground(particles, BackgroundEstimateDefinition);
// 
//   
//   class EventBackgroundWorker{
//     
//     ...
//   };
// 
//   
//   BackgroundEstimate estimate() const;
//   BackgroundEstimate estimate(jet) const;
//   
//   // do we want this:
//   double rho();
//   double rho(jet);
//   //?
//    
// 
//   mutable BackgroundEstimate _estimate;
// 
// 
//   SharedPtr<EventBackgroundWorker> _event_background_worker;
// }
// 
// class BackgroundEstimate{
// 
//   double rho();
//   
//   SharedPtr<EventBackgroundWorker> _event_background_worker;
// 
// private:
//   _rho;
//   _sigma;
//   _rho_m;
//   _sigma_m;
//   _has_rho_m;
//   _has_sigma;
// 
//   // info specific to JMBGE: _reference_jet, mean_area, n_jets_used, n_empty_jets, _empty_area
//   //                         all need to go in the estimate in general
//   // info specific to GMBGE: RectangularGrid, _mean_area (can go either in the the def, the eventBG or the estimate
//   
// };
