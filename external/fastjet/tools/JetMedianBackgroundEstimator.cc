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

#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceStructure.hh"
#include <iostream>
#include <sstream>

using namespace std;

FASTJET_BEGIN_NAMESPACE     // defined in fastjet/internal/base.hh

double BackgroundJetScalarPtDensity::result(const PseudoJet & jet) const {
  // do not include the ghosts in the list of constituents to have a
  // correct behaviour when _pt_power is <= 0
  std::vector<PseudoJet> constituents = (!SelectorIsPureGhost())(jet.constituents());
  double scalar_pt = 0;
  for (unsigned i = 0; i < constituents.size(); i++) {
    scalar_pt += pow(constituents[i].perp(), _pt_power);
  }
  return scalar_pt / jet.area();
}


std::string BackgroundJetScalarPtDensity::description() const {
  ostringstream oss;
  oss << "BackgroundScalarJetPtDensity";
  if (_pt_power != 1.0) oss << " with pt_power = " << _pt_power;
  return oss.str();
}


//----------------------------------------------------------------------
double BackgroundRescalingYPolynomial::result(const PseudoJet & jet) const {
  double y = jet.rap();
  double y2 = y*y;
  double rescaling = _a0 + _a1*y + _a2*y2 + _a3*y2*y + _a4*y2*y2;
  return rescaling;
}

/// allow for warnings
LimitedWarning JetMedianBackgroundEstimator::_warnings;
LimitedWarning JetMedianBackgroundEstimator::_warnings_zero_area;
LimitedWarning JetMedianBackgroundEstimator::_warnings_preliminary;


//---------------------------------------------------------------------
// class JetMedianBackgroundEstimator
// Class to estimate the density of the background per unit area
//---------------------------------------------------------------------

//----------------------------------------------------------------------
// ctors and dtors
//----------------------------------------------------------------------
// ctor that allows to set only the particles later on
JetMedianBackgroundEstimator::JetMedianBackgroundEstimator(const Selector &rho_range,
                                         const JetDefinition &jet_def,
                                         const AreaDefinition &area_def)
  : _rho_range(rho_range), _jet_def(jet_def), _area_def(area_def){

  // initialise things decently
  reset();

  // make a few checks
  _check_jet_alg_good_for_median();
}


//----------------------------------------------------------------------
// ctor from a cluster sequence
//  - csa        the ClusterSequenceArea to use
//  - rho_range  the Selector specifying which jets will be considered
JetMedianBackgroundEstimator::JetMedianBackgroundEstimator( const Selector &rho_range, const ClusterSequenceAreaBase &csa)
  : _rho_range(rho_range), _jet_def(JetDefinition()){

  // initialise things properly
  reset();

  // tell the BGE about the cluster sequence
  set_cluster_sequence(csa);
}


//----------------------------------------------------------------------
// setting a new event
//----------------------------------------------------------------------
// tell the background estimator that it has a new event, composed
// of the specified particles.
void JetMedianBackgroundEstimator::set_particles(const vector<PseudoJet> & particles) {
  // pass an empty seed vector to the full set_particles method to tell it to use
  // default seeds rather than fixed seeds
  vector<int> seed;
  set_particles_with_seed(particles, seed);
}


// tell the background estimator that it has a new event, composed
// of the specified particles and use the supplied seed for the
// generation of ghosts. If the seed is empty, it is ignored.
void JetMedianBackgroundEstimator::set_particles_with_seed(const vector<PseudoJet> & particles, const vector<int> & seed) {
  // make sure that we have been provided a genuine jet definition 
  if (_jet_def.jet_algorithm() == undefined_jet_algorithm)
    throw Error("JetMedianBackgroundEstimator::set_particles can only be called if you set the jet (and area) definition explicitly through the class constructor");

  // initialise things decently (including setting uptodate to false!)
  //reset();

  // cluster the particles
  // 
  // One may argue that it is better to cache the particles and only
  // do the clustering later but clustering the particles now has 2
  // practical advantages:
  //  - it allows us to use only '_included_jets' in all that follows
  //  - it avoids adding another flag to ensure particles are 
  //    clustered only once
  ClusterSequenceArea *csa;
  if (seed.size() == 0) {
    csa = new ClusterSequenceArea(particles, _jet_def, _area_def);
  } else {
    csa = new ClusterSequenceArea(particles, _jet_def, _area_def.with_fixed_seed(seed));
  }
//THREAD-SAFETY-QUESTION: #ifdef FASTJET_HAVE_LIMITED_THREAD_SAFETY
//THREAD-SAFETY-QUESTION:   // before caching thing, lock things down to avoid concurrency issues
//THREAD-SAFETY-QUESTION:   std::lock_guard<std::mutex> guard(_jets_caching_mutex);
//THREAD-SAFETY-QUESTION: #endif
  _included_jets = csa->inclusive_jets();

  // store the CS for later on
  _csi = csa->structure_shared_ptr();
  csa->delete_self_when_unused();

  _set_cache_unavailable();

//THREAD-SAFETY-QUESTION:   // in thread-safe mode, the lock will automatically be released here
}

//----------------------------------------------------------------------
// (re)set the cluster sequence (with area support) to be used by
// future calls to rho() etc. 
//
// \param csa  the cluster sequence area
//
// Pre-conditions: 
//  - one should be able to estimate the "empty area" (i.e. the area
//    not occupied by jets). This is feasible if at least one of the following
//    conditions is satisfied:
//     ( i) the ClusterSequence has explicit ghosts
//     (ii) the range selected has a computable area.
//  - the jet algorithm must be suited for median computation
//    (otherwise a warning will be issues)
//
// Note that selectors with e.g. hardest-jets exclusion do not have
// a well-defined area. For this reasons, it is STRONGLY advised to
// use an area with explicit ghosts.
void JetMedianBackgroundEstimator::set_cluster_sequence(const ClusterSequenceAreaBase & csa) {
  // sanity checks
  //---------------
  //  (i) check that, if there are no explicit ghosts, the selector has a finite area
  if ((!csa.has_explicit_ghosts()) && (!_rho_range.has_finite_area())){
    throw Error("JetMedianBackgroundEstimator: either an area with explicit ghosts (recommended) or a Selector with finite area is needed (to allow for the computation of the empty area)");
  }

//THREAD-SAFETY-QUESTION: #ifdef FASTJET_HAVE_LIMITED_THREAD_SAFETY
//THREAD-SAFETY-QUESTION:   // before caching thing, lock things down to avoid concurrency issues
//THREAD-SAFETY-QUESTION:   std::lock_guard<std::mutex> guard(_jets_caching_mutex);
//THREAD-SAFETY-QUESTION: #endif
  _csi = csa.structure_shared_ptr();

  //  (ii) check the alg is appropriate
  _check_jet_alg_good_for_median();


  // get the initial list of jets
  _included_jets = csa.inclusive_jets();

  _set_cache_unavailable();
  
//THREAD-SAFETY-QUESTION:   // in thread-safe mode, the lock will automatically be released here
}


//----------------------------------------------------------------------
// (re)set the jets (which must have area support) to be used by future
// calls to rho() etc.; for the conditions that must be satisfied
// by the jets, see the Constructor that takes jets.
void JetMedianBackgroundEstimator::set_jets(const vector<PseudoJet> &jets) {
  
  if (! jets.size())
    throw Error("JetMedianBackgroundEstimator::JetMedianBackgroundEstimator: At least one jet is needed to compute the background properties");

  // sanity checks
  //---------------
  //  (o) check that there is an underlying CS shared by all the jets
  if (! (jets[0].has_associated_cluster_sequence()) && (jets[0].has_area()))
    throw Error("JetMedianBackgroundEstimator::JetMedianBackgroundEstimator: the jets used to estimate the background properties must be associated with a valid ClusterSequenceAreaBase");

  SharedPtr<PseudoJetStructureBase> csi_shared = jets[0].structure_shared_ptr();
  ClusterSequenceStructure * csi = dynamic_cast<ClusterSequenceStructure*>(csi_shared.get());
  const ClusterSequenceAreaBase * csab = csi->validated_csab();

  for (unsigned int i=1;i<jets.size(); i++){
    if (! jets[i].has_associated_cluster_sequence()) // area automatic if the next test succeeds
      throw Error("JetMedianBackgroundEstimator::set_jets(...): the jets used to estimate the background properties must be associated with a valid ClusterSequenceAreaBase");

    if (jets[i].structure_shared_ptr().get() != csi_shared.get())
      throw Error("JetMedianBackgroundEstimator::set_jets(...): all the jets used to estimate the background properties must share the same ClusterSequence");
  }

  //  (i) check that, if there are no explicit ghosts, the selector has a finite area
  if ((!csab->has_explicit_ghosts()) && (!_rho_range.has_finite_area())){
    throw Error("JetMedianBackgroundEstimator: either an area with explicit ghosts (recommended) or a Selector with finite area is needed (to allow for the computation of the empty area)");
  }

//THREAD-SAFETY-QUESTION: #ifdef FASTJET_HAVE_LIMITED_THREAD_SAFETY
//THREAD-SAFETY-QUESTION:   // before caching thing, lock things down to avoid concurrency issues
//THREAD-SAFETY-QUESTION:   std::lock_guard<std::mutex> guard(_jets_caching_mutex);
//THREAD-SAFETY-QUESTION: #endif
  _csi = csi_shared;

  //  (ii) check the alg is appropriate
  _check_jet_alg_good_for_median();

  // get the initial list of jets
  _included_jets = jets;

  // ensure recalculation of quantities that need it
  _set_cache_unavailable();
  
//THREAD-SAFETY-QUESTION:   // in thread-safe mode, the lock will automatically be released here
}

//----------------------------------------------------------------------
// retrieving fundamental information
//----------------------------------------------------------------------

// get the full set of background properties
//
// For background estimators using a local ranges, this throws an
//   error (use estimate(jet) instead)
// In the presence of a rescaling, the rescaling factor is not taken
// into account
BackgroundEstimate JetMedianBackgroundEstimator::estimate() const{
  if (_rho_range.takes_reference())
    throw Error("The background estimation is obtained from a selector that takes a reference jet. estimate(PseudoJet) should be used in that case");

  if (!_cache_available) _compute_and_cache_no_overwrite();
  return _cached_estimate;
}

// get the full set of background properties for a given reference jet
// This does not affect the cache
BackgroundEstimate JetMedianBackgroundEstimator::estimate(const PseudoJet &jet) const{
  // first compute an optional rescaling factor
  double rescaling_factor = (_rescaling_class != 0)
    ? (*_rescaling_class)(jet) : 1.0;
  BackgroundEstimate local_estimate;
  
  // adopt a different strategy for ranges taking a reference and others
  if (_rho_range.takes_reference()){
    // we compute the background and rescale it (no caching)
    local_estimate = _compute(jet);
  } else {
    // otherwise, we're in a situation where things can be cached once
    // and for all and then the cache can be used freely
    if (!_cache_available) _compute_and_cache_no_overwrite();
    local_estimate = _cached_estimate;
  }  
  local_estimate.apply_rescaling_factor(rescaling_factor);
  return local_estimate;
}


//------
// get rho, the median background density per unit area
double JetMedianBackgroundEstimator::rho() const {
  if (_rho_range.takes_reference())
    throw Error("The background estimation is obtained from a selector that takes a reference jet. rho(PseudoJet) should be used in that case");

  // we are in a situation where the cache only needs to be computed
  // once, but once it has been computed, we can use it freely.
  if (!_cache_available) _compute_and_cache_no_overwrite();
  return _cached_estimate.rho();
}

// get sigma, the background fluctuations per unit area
double JetMedianBackgroundEstimator::sigma() const {
  if (_rho_range.takes_reference())
    throw Error("The background estimation is obtained from a selector that takes a reference jet. sigma(PseudoJet) should be used in that case");
  if (!_cache_available) _compute_and_cache_no_overwrite();
  return _cached_estimate.sigma();
}

// get rho, the median background density per unit area, locally at
// the position of a given jet.
//
// If the Selector associated with the range takes a reference jet
// (i.e. is relocatable), then for subsequent operations the
// Selector has that jet set as its reference.
double JetMedianBackgroundEstimator::rho(const PseudoJet & jet) {
  // first compute an optional rescaling factor
  double rescaling_factor = (_rescaling_class != 0)
    ? (*_rescaling_class)(jet) : 1.0;
  
  // adopt a different strategy for ranges taking a reference and others
  if (_rho_range.takes_reference()){
    // we compute the background and use it
    BackgroundEstimate estimate = _compute_and_cache_if_needed(jet);
    return rescaling_factor * estimate.rho();
  }

  // otherwise, we're in a situation where things can be cached once
  // and for all and then the cache can be used frely
  if (!_cache_available) _compute_and_cache_no_overwrite();
  return rescaling_factor * _cached_estimate.rho();
}

// get sigma, the background fluctuations per unit area,
// locally at the position of a given jet.
//
// If the Selector associated with the range takes a reference jet
// (i.e. is relocatable), then for subsequent operations the
// Selector has that jet set as its reference.
double JetMedianBackgroundEstimator::sigma(const PseudoJet &jet) {
  // first compute an optional rescaling factor
  double rescaling_factor = (_rescaling_class != 0)
    ? (*_rescaling_class)(jet) : 1.0;
  
  // see "rho(jet)" for a descrtiption of the strategy
  if (_rho_range.takes_reference()){
    BackgroundEstimate estimate = _compute_and_cache_if_needed(jet);
    return rescaling_factor * estimate.sigma();
  }
  // otherwise, cache things once and for all
  if (!_cache_available) _compute_and_cache_no_overwrite();
  return rescaling_factor * _cached_estimate.sigma();
}
 
 
//----------------------------------------------------------------------
// returns rho_m (particle-masses contribution to the 4-vector density)
double JetMedianBackgroundEstimator::rho_m() const {
  if (! has_rho_m()){
    throw Error("JetMediamBackgroundEstimator: rho_m requested but rho_m calculation is disabled (either eplicitly or due to the presence of a jet density class).");
  }
  if (_rho_range.takes_reference()){
    throw Error("The background estimation is obtained from a selector that takes a reference jet. rho_m(PseudoJet) should be used in that case");
  }
  if (!_cache_available) _compute_and_cache_no_overwrite();
  return _cached_estimate.rho_m();
}


//----------------------------------------------------------------------
// returns sigma_m (particle-masses contribution to the 4-vector
// density); must be multipled by sqrt(area) to get fluctuations
// for a region of a given area.
double JetMedianBackgroundEstimator::sigma_m() const{
  if (! has_rho_m()){
    throw Error("JetMediamBackgroundEstimator: sigma_m requested but rho_m/sigma_m calculation is disabled (either explicitly or due to the presence of a jet density class).");
  }
  if (_rho_range.takes_reference())
    throw Error("The background estimation is obtained from a selector that takes a reference jet. sigma_m(PseudoJet) should be used in that case");
  if (!_cache_available) _compute_and_cache_no_overwrite();
  return _cached_estimate.sigma_m();
}

//----------------------------------------------------------------------
// returns rho_m locally at the position of a given jet. As for
// rho(jet), it is non-const.
double JetMedianBackgroundEstimator::rho_m(const PseudoJet & jet)  {
  // first compute an optional rescaling factor
  double rescaling_factor = (_rescaling_class != 0)
    ? (*_rescaling_class)(jet) : 1.0;
  
  // see "rho(jet)" for a descrtiption of the strategy
  if (_rho_range.takes_reference()){
    BackgroundEstimate estimate = _compute_and_cache_if_needed(jet);
    return rescaling_factor * estimate.rho_m();
  }  
  // otherwise, cache things once and for all
  if (!_cache_available) _compute_and_cache_no_overwrite();
  return rescaling_factor * _cached_estimate.rho_m();
}


//----------------------------------------------------------------------
// returns sigma_m locally at the position of a given jet. As for
// rho(jet), it is non-const.
double JetMedianBackgroundEstimator::sigma_m(const PseudoJet & jet){
  // first compute an optional rescaling factor
  double rescaling_factor = (_rescaling_class != 0)
    ? (*_rescaling_class)(jet) : 1.0;
  
  // see "rho(jet)" for a descrtiption of the strategy
  if (_rho_range.takes_reference()){
    BackgroundEstimate estimate = _compute_and_cache_if_needed(jet);
    return rescaling_factor * estimate.sigma_m();
  }
  // otherwise, cache things once and for all
  if (!_cache_available) _compute_and_cache_no_overwrite();
  return rescaling_factor * _cached_estimate.sigma_m();
}


//----------------------------------------------------------------
/// Returns the mean area of the jets used to actually compute the
/// background properties in the last call of rho() or sigma()
/// If the configuration has changed in the meantime, throw an error.
double JetMedianBackgroundEstimator::mean_area() const{
  // if the selector takes a reference, we need to use the cache 
  if (_rho_range.takes_reference()){
    // lock to make sure no other thread interferes w the caching
    _lock_if_needed();
    if (!_cache_available){
      _unlock_if_needed();
      throw Error("Calls to JetMedianBackgroundEstimator::mean_area() in cases where the background estimation uses a selector that takes a reference jet need to call a method that fills the cached estimate (rho(jet), sigma(jet), ...).");
    }
    double return_value = _cached_estimate.mean_area();
    _unlock_if_needed();
    return return_value;
  }
  if (!_cache_available) _compute_and_cache_no_overwrite();
  return _cached_estimate.mean_area();
}

/// returns the number of jets used to actually compute the
/// background properties in the last call of rho() or sigma()
/// If the configuration has changed in the meantime, throw an error.
unsigned int JetMedianBackgroundEstimator::n_jets_used() const{
  // if the selector takes a reference, we need to use the cache 
  if (_rho_range.takes_reference()){
    // lock to make sure no other thread interferes w the caching
    _lock_if_needed();
    if (!_cache_available){
      _unlock_if_needed();
      throw Error("Calls to JetMedianBackgroundEstimator::n_jets_used() in cases where the background estimation uses a selector that takes a reference jet need to call a method that fills the cached estimate (rho(jet), sigma(jet), ...).");
    }
    unsigned int return_value = _cached_estimate.extras<JetMedianBackgroundEstimator>().n_jets_used();
    _unlock_if_needed();
    return return_value;
  }
  if (!_cache_available) _compute_and_cache_no_overwrite();
  return _cached_estimate.extras<JetMedianBackgroundEstimator>().n_jets_used();
}

/// returns the jets used to actually compute the background
/// properties
std::vector<PseudoJet> JetMedianBackgroundEstimator::jets_used() const{
  vector<PseudoJet> tmp_jets;
  
  // if the selector takes a reference, we need to use the cache 
  if (_rho_range.takes_reference()){
    // lock to make sure no other thread interferes w the caching
    _lock_if_needed();
    if (!_cache_available){
      _unlock_if_needed();
      throw Error("Calls to JetMedianBackgroundEstimator::jets_used() in cases where the background estimation uses a selector that takes a reference jet need to call a method that fills the cached estimate (rho(jet), sigma(jet), ...).");
    }
    PseudoJet reference_jet = _cached_estimate.extras<JetMedianBackgroundEstimator>().reference_jet();
    _unlock_if_needed();
    Selector local_rho_range = _rho_range;
    local_rho_range.set_reference(reference_jet);
    tmp_jets = local_rho_range(_included_jets);
  } else {
    if (!_cache_available) _compute_and_cache_no_overwrite();
    tmp_jets = _rho_range(_included_jets);
  }
  
  std::vector<PseudoJet> used_jets;
  for (unsigned int i=0; i<tmp_jets.size(); i++){
    if (tmp_jets[i].area()>0) used_jets.push_back(tmp_jets[i]);
  }
  return used_jets;
}

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
double JetMedianBackgroundEstimator::empty_area() const{
  // if the selector takes a reference, we need to use the cache 
  if (_rho_range.takes_reference()){
    // lock to make sure no other thread interferes w the caching
    _lock_if_needed();
    if (!_cache_available){
      _unlock_if_needed();
      throw Error("Calls to JetMedianBackgroundEstimator::empty_area() in cases where the background estimation uses a selector that takes a reference jet need to call a method that fills the cached estimate (rho(jet), sigma(jet), ...).");
    }
    double return_value = _cached_estimate.extras<JetMedianBackgroundEstimator>().empty_area();
    _unlock_if_needed();
    return return_value;
  }
  if (!_cache_available) _compute_and_cache_no_overwrite();
  return _cached_estimate.extras<JetMedianBackgroundEstimator>().empty_area();
}

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
double JetMedianBackgroundEstimator::n_empty_jets() const{
  // if the selector takes a reference, we need to use the cache 
  if (_rho_range.takes_reference()){
    // lock to make sure no other thread interferes w the caching
    _lock_if_needed();
    if (!_cache_available){
      _unlock_if_needed();
      throw Error("Calls to JetMedianBackgroundEstimator::n_empty_jets() in cases where the background estimation uses a selector that takes a reference jet need to call a method that fills the cached estimate (rho(jet), sigma(jet), ...).");
    }
    double return_value = _cached_estimate.extras<JetMedianBackgroundEstimator>().n_empty_jets();
    _unlock_if_needed();
    return return_value;
  }
  if (!_cache_available) _compute_and_cache_no_overwrite();
  return _cached_estimate.extras<JetMedianBackgroundEstimator>().n_empty_jets();
}
 

//----------------------------------------------------------------------
// configuring behaviour
//----------------------------------------------------------------------
// reset to default values
// 
// set the variou options to their default values
void JetMedianBackgroundEstimator::reset(){
  // set the remaining default parameters
  set_use_area_4vector();  // true by default
  set_provide_fj2_sigma(false);

  _enable_rho_m = true;

  // reset the computed values
  _included_jets.clear();

  _jet_density_class = 0; // null pointer
  _rescaling_class = 0;   // null pointer

  _set_cache_unavailable();
}


// Set a pointer to a class that calculates the quantity whose
// median will be calculated; if the pointer is null then pt/area
// is used (as occurs also if this function is not called).
void JetMedianBackgroundEstimator::set_jet_density_class(const FunctionOfPseudoJet<double> * jet_density_class_in) {
  _warnings_preliminary.warn("JetMedianBackgroundEstimator::set_jet_density_class: density classes are still preliminary in FastJet 3.1. Their interface may differ in future releases (without guaranteeing backward compatibility). Note that since FastJet 3.1, rho_m and sigma_m are accessible direclty in JetMedianBackgroundEstimator and GridMedianBackgroundEstimator(with no need for a density class).");
  _jet_density_class = jet_density_class_in;
  _set_cache_unavailable();
}



//----------------------------------------------------------------------
// description
//----------------------------------------------------------------------
string JetMedianBackgroundEstimator::description() const { 
  ostringstream desc;
  desc << "JetMedianBackgroundEstimator, using " << _jet_def.description() 
       << " with " << _area_def.description() 
       << " and selecting jets with " << _rho_range.description();
  return desc.str();
}       



//----------------------------------------------------------------------
// computation of the background properties
//----------------------------------------------------------------------
// do the actual job
BackgroundEstimate JetMedianBackgroundEstimator::_compute(const PseudoJet &jet) const {
  // prepare a local structure to hold temporarily the results
  // (by design, this comes with default values of 0 for each property)
  BackgroundEstimate local_estimate;

  // check if the clustersequence is still valid
  _check_csa_alive();

  local_estimate.set_has_sigma(has_sigma());
  local_estimate.set_has_rho_m(has_rho_m());

  // structure to hold the extra info associated w this BGE (the call
  // below initialises everything to 0)
  Extras * extras = new Extras;
  local_estimate.set_extras(extras);
  extras->set_reference_jet(jet);

  // fill the vector of pt/area (or the quantity from the jet density class) 
  //  - in the range
  vector<double> vector_for_median_pt;
  vector<double> vector_for_median_dt;
  double total_area  = 0.0;
  
  // apply the selector to the included jets
  vector<PseudoJet> selected_jets;
  if (_rho_range.takes_reference()){
    Selector local_rho_range = _rho_range;
    selected_jets = local_rho_range.set_reference(jet)(_included_jets);
  } else {
    selected_jets = _rho_range(_included_jets);
  }
  
  // compute the pt/area for the selected jets
  double median_input_pt, median_input_dt=0.0;
  BackgroundJetPtMDensity m_density;
  bool do_rho_m = has_rho_m();
  unsigned int njets_used = 0;
  for (unsigned i = 0; i < selected_jets.size(); i++) {
    const PseudoJet & current_jet = selected_jets[i];

    double this_area = (_use_area_4vector) ? current_jet.area_4vector().perp() : current_jet.area(); 
    if (this_area>0){
      // for the pt component, we either use pt or the user-provided
      // density class
      if (_jet_density_class == 0) {
        median_input_pt = current_jet.perp()/this_area;
      } else {
        median_input_pt = (*_jet_density_class)(current_jet);
      }

      // handle the rho_m part if requested
      // note that we're using the scalar area as a normalisation inside the
      // density class!
      if (do_rho_m) 
        median_input_dt = m_density(current_jet);
    
      // perform rescaling if needed
      if (_rescaling_class != 0) {
        double resc = (*_rescaling_class)(current_jet);
        median_input_pt /= resc;
        median_input_dt /= resc;
      }
      
      // store the result for future computation of the median
      vector_for_median_pt.push_back(median_input_pt);
      if (do_rho_m) 
        vector_for_median_dt.push_back(median_input_dt);

      total_area  += this_area;
      njets_used++;
    } else {
      _warnings_zero_area.warn("JetMedianBackgroundEstimator::_compute(...): discarded jet with zero area. Zero-area jets may be due to (i) too large a ghost area (ii) a jet being outside the ghost range (iii) the computation not being done using an appropriate algorithm (kt;C/A).");
    }
  }
  
  // there is nothing inside our region, so answer will always be zero
  if (vector_for_median_pt.size() == 0) {
    // record that the computation has been performed  
    return local_estimate;
  }

  // determine the number of empty jets
  // If we have explicit ghosts, this is 0 (i.e. the default)
  const ClusterSequenceAreaBase * csab = (dynamic_cast<ClusterSequenceStructure*>(_csi.get()))->validated_csab();
  if (! (csab->has_explicit_ghosts())) {
    if (_rho_range.takes_reference()){
      Selector local_rho_range = _rho_range;
      local_rho_range.set_reference(jet);
      extras->set_empty_area  (csab->empty_area(local_rho_range));
      extras->set_n_empty_jets(csab->n_empty_jets(local_rho_range));
    } else {
      extras->set_empty_area  (csab->empty_area(_rho_range));
      extras->set_n_empty_jets(csab->n_empty_jets(_rho_range));
    }
  }

  extras->set_n_jets_used(njets_used);
  double total_njets = extras->n_jets_used() + extras->n_empty_jets();
  total_area  += extras->empty_area();

  double rho_tmp, stand_dev;
  _median_and_stddev(vector_for_median_pt, extras->n_empty_jets(),
                     rho_tmp, stand_dev, _provide_fj2_sigma);
  local_estimate.set_rho(rho_tmp);
  
  // process and store the results (_rho was already stored above)
  local_estimate.set_mean_area(total_area / total_njets);
  local_estimate.set_sigma(stand_dev * sqrt( max(0.0, local_estimate.mean_area()) ));

  // compute the rho_m part now
  if (do_rho_m){
    _median_and_stddev(vector_for_median_dt, extras->n_empty_jets(),
                       rho_tmp, stand_dev, 
		       _provide_fj2_sigma);
    local_estimate.set_rho_m(rho_tmp);
    local_estimate.set_sigma_m(stand_dev * sqrt( max(0.0, local_estimate.mean_area()) ));
  }

  return local_estimate;
}


void JetMedianBackgroundEstimator::_cache_no_overwrite(const BackgroundEstimate &estimate) const {
  /// this is meant to be called if the selector is not local
  assert(!(_rho_range.takes_reference()));
  
  // we need to write to the cache, so set a lock if needed
  _lock_if_needed();

  // we only need to write if someone else did not do it earlier
  //
  // Doing this avoids potential "read" problems when two threads try
  // to compute the cache at the same time. The first one may return
  // and try to read the result when the second actually writes. The
  // lines below guarantee that the first thread would have set
  // _cache_available to true before releasing the lock and therefore
  // the second thread will not attempt to write
  if (!_cache_available){
    _cached_estimate = estimate;
    _cache_available = true;
  }

  // release the lock
  _unlock_if_needed();
}

void JetMedianBackgroundEstimator::_compute_and_cache_no_overwrite() const {
  /// this is meant to be called if the selector is not local
  assert(!(_rho_range.takes_reference()));
  
  // get the result (for a dummy PseudoJet which will anyway not be used)
  // and cache it
  _cache_no_overwrite(_compute(PseudoJet()));
}
 
void JetMedianBackgroundEstimator::_cache(const BackgroundEstimate &estimate) const {
  /// this is meant to be called if the selector is local
  assert(_rho_range.takes_reference());

  // we need to write to the cache, so set a lock if needed
  _lock_if_needed();

  // if we overwrite the cache, we need to make sure that other
  // parts of the code use r/w accesses that respect the lock that
  // we have acquired.
  //
  // In practice, this will only be used in the case wheere we have
  // a local selector, in queries that require access to the cache.
  _cached_estimate = estimate;
  _cache_available = true;    

  // release the lock
  _unlock_if_needed();
}

BackgroundEstimate JetMedianBackgroundEstimator::_compute_and_cache_if_needed(const PseudoJet &jet) const {
  /// this is meant to be called if the selector is local
  assert(_rho_range.takes_reference());

  BackgroundEstimate local_estimate;
  
  _lock_if_needed();
  if ((_cache_available) && (_cached_estimate.extras<JetMedianBackgroundEstimator>().reference_jet() == jet)){
    local_estimate = _cached_estimate;
    _unlock_if_needed();
    return local_estimate;
  }
  _unlock_if_needed();

  local_estimate = _compute(jet);
  _cache(local_estimate);
  return local_estimate;
}

// check that the underlying structure is still alive;
// throw an error otherwise
void JetMedianBackgroundEstimator::_check_csa_alive() const{
  ClusterSequenceStructure* csa = dynamic_cast<ClusterSequenceStructure*>(_csi.get());
  if (csa == 0) {
    throw Error("JetMedianBackgroundEstimator: there is no cluster sequence associated with the JetMedianBackgroundEstimator");
  }
  if (! csa->has_associated_cluster_sequence())
    throw Error("JetMedianBackgroundEstimator: modifications are no longer possible as the underlying ClusterSequence has gone out of scope");
}


// check that the algorithm used for the clustering is suitable for
// background estimation (i.e. either kt or C/A).
// Issue a warning otherwise
void JetMedianBackgroundEstimator::_check_jet_alg_good_for_median() const{
  const JetDefinition * jet_def = &_jet_def;

  // if no explicit jet def has been provided, fall back on the
  // cluster sequence
  if (_jet_def.jet_algorithm() == undefined_jet_algorithm){
    const ClusterSequence * cs = dynamic_cast<ClusterSequenceStructure*>(_csi.get())->validated_cs();
    jet_def = &(cs->jet_def());
  }

  if (jet_def->jet_algorithm() != kt_algorithm
      && jet_def->jet_algorithm() != cambridge_algorithm
      && jet_def->jet_algorithm() != cambridge_for_passive_algorithm) {
    _warnings.warn("JetMedianBackgroundEstimator: jet_def being used may not be suitable for estimating diffuse backgrounds (good alternatives are kt, cam)");
  }
}


FASTJET_END_NAMESPACE


