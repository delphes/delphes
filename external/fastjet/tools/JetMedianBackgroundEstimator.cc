//FJSTARTHEADER
// $Id: JetMedianBackgroundEstimator.cc 3517 2014-08-01 14:23:13Z soyez $
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

#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/ClusterSequenceStructure.hh>
#include <iostream>
#include <sstream>

FASTJET_BEGIN_NAMESPACE     // defined in fastjet/internal/base.hh

using namespace std;

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
  // make sure that we have been provided a genuine jet definition 
  if (_jet_def.jet_algorithm() == undefined_jet_algorithm)
    throw Error("JetMedianBackgroundEstimator::set_particles can only be called if you set the jet (and area) definition explicitly through the class constructor");

  // initialise things decently (including setting uptodate to false!)
  //reset();
  _uptodate=false;

  // cluster the particles
  // 
  // One may argue that it is better to cache the particles and only
  // do the clustering later but clustering the particles now has 2
  // practical advantages:
  //  - it allows to une only '_included_jets' in all that follows
  //  - it avoids adding another flag to ensure particles are 
  //    clustered only once
  ClusterSequenceArea *csa = new ClusterSequenceArea(particles, _jet_def, _area_def);
  _included_jets = csa->inclusive_jets();

  // store the CS for later on
  _csi = csa->structure_shared_ptr();
  csa->delete_self_when_unused();
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
  _csi = csa.structure_shared_ptr();

  // sanity checks
  //---------------
  //  (i) check the alg is appropriate
  _check_jet_alg_good_for_median();

  //  (ii) check that, if there are no explicit ghosts, the selector has a finite area
  if ((!csa.has_explicit_ghosts()) && (!_rho_range.has_finite_area())){
    throw Error("JetMedianBackgroundEstimator: either an area with explicit ghosts (recommended) or a Selector with finite area is needed (to allow for the computation of the empty area)");
  }

  // get the initial list of jets
  _included_jets = csa.inclusive_jets();

  _uptodate = false;
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

  _csi = jets[0].structure_shared_ptr();
  ClusterSequenceStructure * csi = dynamic_cast<ClusterSequenceStructure*>(_csi());
  const ClusterSequenceAreaBase * csab = csi->validated_csab();

  for (unsigned int i=1;i<jets.size(); i++){
    if (! jets[i].has_associated_cluster_sequence()) // area automatic if the next test succeeds
      throw Error("JetMedianBackgroundEstimator::set_jets(...): the jets used to estimate the background properties must be associated with a valid ClusterSequenceAreaBase");

    if (jets[i].structure_shared_ptr().get() != _csi.get())
      throw Error("JetMedianBackgroundEstimator::set_jets(...): all the jets used to estimate the background properties must share the same ClusterSequence");
  }

  //  (i) check the alg is appropriate
  _check_jet_alg_good_for_median();

  //  (ii) check that, if there are no explicit ghosts, the selector has a finite area
  if ((!csab->has_explicit_ghosts()) && (!_rho_range.has_finite_area())){
    throw Error("JetMedianBackgroundEstimator: either an area with explicit ghosts (recommended) or a Selector with finite area is needed (to allow for the computation of the empty area)");
  }


  // get the initial list of jets
  _included_jets = jets;

  // ensure recalculation of quantities that need it
  _uptodate = false;
}


//----------------------------------------------------------------------
// retrieving fundamental information
//----------------------------------------------------------------

// get rho, the median background density per unit area
double JetMedianBackgroundEstimator::rho() const {
  if (_rho_range.takes_reference())
    throw Error("The background estimation is obtained from a selector that takes a reference jet. rho(PseudoJet) should be used in that case");
  _recompute_if_needed();
  return _rho;
}

// get sigma, the background fluctuations per unit area
double JetMedianBackgroundEstimator::sigma() const {
  if (_rho_range.takes_reference())
    throw Error("The background estimation is obtained from a selector that takes a reference jet. rho(PseudoJet) should be used in that case");
  _recompute_if_needed();
  return _sigma;
}

// get rho, the median background density per unit area, locally at
// the position of a given jet.
//
// If the Selector associated with the range takes a reference jet
// (i.e. is relocatable), then for subsequent operations the
// Selector has that jet set as its reference.
double JetMedianBackgroundEstimator::rho(const PseudoJet & jet) {
  _recompute_if_needed(jet);
  double our_rho = _rho;
  if (_rescaling_class != 0) { 
    our_rho *= (*_rescaling_class)(jet);
  }
  return our_rho;
}

// get sigma, the background fluctuations per unit area,
// locally at the position of a given jet.
//
// If the Selector associated with the range takes a reference jet
// (i.e. is relocatable), then for subsequent operations the
// Selector has that jet set as its reference.
double JetMedianBackgroundEstimator::sigma(const PseudoJet &jet) {
  _recompute_if_needed(jet);
  double our_sigma = _sigma;
  if (_rescaling_class != 0) { 
    our_sigma *= (*_rescaling_class)(jet);
  }
  return our_sigma;
}


//----------------------------------------------------------------------
// returns rho_m (particle-masses contribution to the 4-vector density)
double JetMedianBackgroundEstimator::rho_m() const {
  if (! has_rho_m()){
    throw Error("JetMediamBackgroundEstimator: rho_m requested but rho_m calculation is disabled (either eplicitly or due to the presence of a jet density class).");
  }
  if (_rho_range.takes_reference())
    throw Error("The background estimation is obtained from a selector that takes a reference jet. rho(PseudoJet) should be used in that case");
  _recompute_if_needed();
  return _rho_m;
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
    throw Error("The background estimation is obtained from a selector that takes a reference jet. rho(PseudoJet) should be used in that case");
  _recompute_if_needed();
  return _sigma_m;
}

//----------------------------------------------------------------------
// returns rho_m locally at the position of a given jet. As for
// rho(jet), it is non-const.
double JetMedianBackgroundEstimator::rho_m(const PseudoJet & jet)  {
  _recompute_if_needed(jet);
  double our_rho = _rho_m;
  if (_rescaling_class != 0) { 
    our_rho *= (*_rescaling_class)(jet);
  }
  return our_rho;
}


//----------------------------------------------------------------------
// returns sigma_m locally at the position of a given jet. As for
// rho(jet), it is non-const.
double JetMedianBackgroundEstimator::sigma_m(const PseudoJet & jet){
  _recompute_if_needed(jet);
  double our_sigma = _sigma_m;
  if (_rescaling_class != 0) { 
    our_sigma *= (*_rescaling_class)(jet);
  }
  return our_sigma;
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
  _rho = _sigma = 0.0;
  _rho_m = _sigma_m = 0.0;
  _n_jets_used = _n_empty_jets = 0;
  _empty_area = _mean_area = 0.0;

  _included_jets.clear();

  _jet_density_class = 0; // null pointer
  _rescaling_class = 0;   // null pointer

  _uptodate = false;
}


// Set a pointer to a class that calculates the quantity whose
// median will be calculated; if the pointer is null then pt/area
// is used (as occurs also if this function is not called).
void JetMedianBackgroundEstimator::set_jet_density_class(const FunctionOfPseudoJet<double> * jet_density_class_in) {
  _warnings_preliminary.warn("JetMedianBackgroundEstimator::set_jet_density_class: density classes are still preliminary in FastJet 3.1. Their interface may differ in future releases (without guaranteeing backward compatibility). Note that since FastJet 3.1, rho_m and sigma_m are accessible direclty in JetMedianBackgroundEstimator and GridMedianBackgroundEstimator(with no need for a density class).");
  _jet_density_class = jet_density_class_in;
  _uptodate = false;
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
// for estimation using a relocatable selector (i.e. local range)
// this allows to set its position. Note that this HAS to be called
// before any attempt to compute the background properties
void JetMedianBackgroundEstimator::_recompute_if_needed(const PseudoJet &jet){
  // if the range is relocatable, handles its relocation
  if (_rho_range.takes_reference()){
    // check that the reference is not the same as the previous one
    // (would avoid an unnecessary recomputation)
    if (jet == _current_reference) return;

    // relocate the range and make sure things get recomputed the next
    // time one tries to get some information
    _rho_range.set_reference(jet);
    _uptodate=false;
  }

  _recompute_if_needed();
}

// do the actual job
void JetMedianBackgroundEstimator::_compute() const {
  // check if the clustersequence is still valid
  _check_csa_alive();

  // fill the vector of pt/area (or the quantity from the jet density class) 
  //  - in the range
  vector<double> vector_for_median_pt;
  vector<double> vector_for_median_dt;
  double total_area  = 0.0;
  _n_jets_used = 0;

  // apply the selector to the included jets
  vector<PseudoJet> selected_jets = _rho_range(_included_jets);

  // compute the pt/area for the selected jets
  double median_input_pt, median_input_dt=0.0;
  BackgroundJetPtMDensity m_density;
  bool do_rho_m = has_rho_m();
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
	double resc = (*_rescaling_class)(current_jet);;
        median_input_pt /= resc;
        median_input_dt /= resc;
      }
      
      // store the result for future computation of the median
      vector_for_median_pt.push_back(median_input_pt);
      if (do_rho_m) 
	vector_for_median_dt.push_back(median_input_dt);

      total_area  += this_area;
      _n_jets_used++;
    } else {
      _warnings_zero_area.warn("JetMedianBackgroundEstimator::_compute(...): discarded jet with zero area. Zero-area jets may be due to (i) too large a ghost area (ii) a jet being outside the ghost range (iii) the computation not being done using an appropriate algorithm (kt;C/A).");
    }
  }
  
  // there is nothing inside our region, so answer will always be zero
  if (vector_for_median_pt.size() == 0) {
    _rho        = 0.0;
    _sigma      = 0.0;
    _rho_m      = 0.0;
    _sigma_m    = 0.0;
    _mean_area  = 0.0;
    return;
  }

  // determine the number of empty jets
  const ClusterSequenceAreaBase * csab = (dynamic_cast<ClusterSequenceStructure*>(_csi()))->validated_csab();
  if (csab->has_explicit_ghosts()) {
    _empty_area = 0.0;
    _n_empty_jets = 0;
  } else {
    _empty_area = csab->empty_area(_rho_range);
    _n_empty_jets = csab->n_empty_jets(_rho_range);
  }

  double total_njets = _n_jets_used + _n_empty_jets;
  total_area  += _empty_area;

  double stand_dev;
  _median_and_stddev(vector_for_median_pt, _n_empty_jets, _rho, stand_dev, 
                     _provide_fj2_sigma);

  // process and store the results (_rho was already stored above)
  _mean_area  = total_area / total_njets;
  _sigma      = stand_dev * sqrt(_mean_area);

  // compute the rho_m part now
  if (do_rho_m){
    _median_and_stddev(vector_for_median_dt, _n_empty_jets, _rho_m, stand_dev, 
		       _provide_fj2_sigma);
    _sigma_m = stand_dev * sqrt(_mean_area);
  }

  // record that the computation has been performed  
  _uptodate = true;
}



// check that the underlying structure is still alive;
// throw an error otherwise
void JetMedianBackgroundEstimator::_check_csa_alive() const{
  ClusterSequenceStructure* csa = dynamic_cast<ClusterSequenceStructure*>(_csi());
  if (csa == 0) {
    throw Error("JetMedianBackgroundEstimator: there is no cluster sequence associated with the JetMedianBackgroundEstimator");
  }
  if (! dynamic_cast<ClusterSequenceStructure*>(_csi())->has_associated_cluster_sequence())
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
    const ClusterSequence * cs = dynamic_cast<ClusterSequenceStructure*>(_csi())->validated_cs();
    jet_def = &(cs->jet_def());
  }

  if (jet_def->jet_algorithm() != kt_algorithm
      && jet_def->jet_algorithm() != cambridge_algorithm
      && jet_def->jet_algorithm() != cambridge_for_passive_algorithm) {
    _warnings.warn("JetMedianBackgroundEstimator: jet_def being used may not be suitable for estimating diffuse backgrounds (good alternatives are kt, cam)");
  }
}



FASTJET_END_NAMESPACE


