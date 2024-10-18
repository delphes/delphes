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


#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


//----------------------------------------------------------------------
// setting a new event
//----------------------------------------------------------------------
// tell the background estimator that it has a new event, composed
// of the specified particles.
void GridMedianBackgroundEstimator::set_particles(const vector<PseudoJet> & particles) {
  vector<double> scalar_pt(n_tiles(), 0.0);

  assert(all_tiles_equal_area());
  //assert(n_good_tiles() == n_tiles()); // not needed now that we have an implementation

  _cached_estimate.reset();
  _cached_estimate.set_has_sigma(true);
  _cached_estimate.set_mean_area(mean_tile_area());
  
  // check if we need to compute only rho or both rho and rho_m
  if (_enable_rho_m){
    // both rho and rho_m
    //
    // this requires a few other variables
    vector<double> scalar_dt(n_tiles(), 0.0);
    double pt, dt;
    for (unsigned i = 0; i < particles.size(); i++) {
      int j = tile_index(particles[i]);
      if (j >= 0){
	pt = particles[i].pt();
	dt = particles[i].mt() - pt;
	if (_rescaling_class == 0){
	  scalar_pt[j] += pt;
	  scalar_dt[j] += dt;
	} else {
	  double r = (*_rescaling_class)(particles[i]);
	  scalar_pt[j] += pt/r;
	  scalar_dt[j] += dt/r;
	}
      }
    }
    // sort things for _percentile
    sort(scalar_dt.begin(), scalar_dt.end());

    // compute rho_m and sigma_m (see comment below for the
    // normaliosation of sigma)
    double p50 = _percentile(scalar_dt, 0.5);
    _cached_estimate.set_has_rho_m(true);
    _cached_estimate.set_rho_m(p50 / mean_tile_area());
    _cached_estimate.set_sigma_m((p50-_percentile(scalar_dt, (1.0-0.6827)/2.0))/sqrt(mean_tile_area()));
  } else {
    // only rho
    //fill(_scalar_pt.begin(), _scalar_pt.end(), 0.0);
    for (unsigned i = 0; i < particles.size(); i++) {
      int j = tile_index(particles[i]);
      if (j >= 0){
	if (_rescaling_class == 0){
	  scalar_pt[j] += particles[i].pt();
	} else {
	  scalar_pt[j] += particles[i].pt()/(*_rescaling_class)(particles[i]);
	}
      }
    }
  }

  // if there are some "bad" tiles, then we need to exclude them from
  // the calculation of the median. We'll do this by condensing the
  // scalar_pt vector down to just the values for the tiles that are
  // good.
  //
  // tested answers look right in "issue" 2014-08-08-testing-rect-grid
  if (n_good_tiles() != n_tiles()) {
    int newn = 0;
    for (unsigned i = 0; i < scalar_pt.size(); i++) {
      if (tile_is_good(i)) {
        // clang gets confused with the SharedPtr swap if we don't
        // have std:: here
        std::swap(scalar_pt[i],scalar_pt[newn]);
        newn++;
      }
    }
    scalar_pt.resize(newn);
  }

  // in all cases, carry on with the computation of rho
  // 
  // first sort
  sort(scalar_pt.begin(), scalar_pt.end());

  // then compute rho
  //
  // watch out: by definition, our sigma is the standard deviation of
  // the pt density multiplied by the square root of the cell area
  double p50 = _percentile(scalar_pt, 0.5);
  _cached_estimate.set_rho(p50 / mean_tile_area());
  _cached_estimate.set_sigma((p50-_percentile(scalar_pt, (1.0-0.6827)/2.0))/sqrt(mean_tile_area()));

  _cache_available = true;
}


//----------------------------------------------------------------------
// retrieving fundamental information
//----------------------------------------------------------------------

// get the full set of background properties
BackgroundEstimate GridMedianBackgroundEstimator::estimate() const{
  verify_particles_set();
  return _cached_estimate;  
}
 
// get the full set of background properties for a given reference jet
BackgroundEstimate GridMedianBackgroundEstimator::estimate(const PseudoJet &jet) const{
  verify_particles_set();
  if (_rescaling_class == 0)
    return _cached_estimate;
  
  BackgroundEstimate local_estimate = _cached_estimate;
  local_estimate.apply_rescaling_factor((*_rescaling_class)(jet));
  return local_estimate;
}


// get rho, the median background density per unit area
double GridMedianBackgroundEstimator::rho() const {
  verify_particles_set();
  return _cached_estimate.rho();
}


//----------------------------------------------------------------------
// get sigma, the background fluctuations per unit area; must be
// multipled by sqrt(area) to get fluctuations for a region of a
// given area.
double GridMedianBackgroundEstimator::sigma() const{
  verify_particles_set();
  return _cached_estimate.sigma(); 
}

//----------------------------------------------------------------------
// get rho, the background density per unit area, locally at the
// position of a given jet. Note that this is not const, because a
// user may then wish to query other aspects of the background that
// could depend on the position of the jet last used for a rho(jet)
// determination.
double GridMedianBackgroundEstimator::rho(const PseudoJet & jet)  {
  //verify_particles_set();
  double rescaling = (_rescaling_class == 0) ? 1.0 : (*_rescaling_class)(jet);
  return rescaling*rho();
}


//----------------------------------------------------------------------
// get sigma, the background fluctuations per unit area, locally at
// the position of a given jet. As for rho(jet), it is non-const.
double GridMedianBackgroundEstimator::sigma(const PseudoJet & jet){
  //verify_particles_set();
  double rescaling = (_rescaling_class == 0) ? 1.0 : (*_rescaling_class)(jet);
  return rescaling*sigma();
}

//----------------------------------------------------------------------
// returns rho_m (particle-masses contribution to the 4-vector density)
double GridMedianBackgroundEstimator::rho_m() const {
  if (! _enable_rho_m){
    throw Error("GridMediamBackgroundEstimator: rho_m requested but rho_m calculation has been disabled.");
  }
  verify_particles_set();
  return _cached_estimate.rho_m();
}


//----------------------------------------------------------------------
// returns sigma_m (particle-masses contribution to the 4-vector
// density); must be multipled by sqrt(area) to get fluctuations
// for a region of a given area.
double GridMedianBackgroundEstimator::sigma_m() const{
  if (! _enable_rho_m){
    throw Error("GridMediamBackgroundEstimator: sigma_m requested but rho_m/sigma_m calculation has been disabled.");
  }
  verify_particles_set();
  return _cached_estimate.sigma_m(); 
}

//----------------------------------------------------------------------
// returns rho_m locally at the position of a given jet. As for
// rho(jet), it is non-const.
double GridMedianBackgroundEstimator::rho_m(const PseudoJet & jet)  {
  //verify_particles_set();
  double rescaling = (_rescaling_class == 0) ? 1.0 : (*_rescaling_class)(jet);
  return rescaling*rho_m();
}


//----------------------------------------------------------------------
// returns sigma_m locally at the position of a given jet. As for
// rho(jet), it is non-const.
double GridMedianBackgroundEstimator::sigma_m(const PseudoJet & jet){
  //verify_particles_set();
  double rescaling = (_rescaling_class == 0) ? 1.0 : (*_rescaling_class)(jet);
  return rescaling*sigma_m();
}

//----------------------------------------------------------------------
// verify that particles have been set and throw an error if not
void GridMedianBackgroundEstimator::verify_particles_set() const {
  if (!_cache_available) throw Error("GridMedianBackgroundEstimator::rho() or sigma() called without particles having been set");
}


//----------------------------------------------------------------------
// description
//----------------------------------------------------------------------
string GridMedianBackgroundEstimator::description() const { 
  ostringstream desc;
  desc << "GridMedianBackgroundEstimator, with " << RectangularGrid::description();
  return desc.str();
}       


//----------------------------------------------------------------------
// configuring the behaviour
//----------------------------------------------------------------------
// Set a pointer to a class that calculates the rescaling factor as
// a function of the jet (position). Note that the rescaling factor
// is used both in the determination of the "global" rho (the pt/A
// of each jet is divided by this factor) and when asking for a
// local rho (the result is multiplied by this factor).
//
// The BackgroundRescalingYPolynomial class can be used to get a
// rescaling that depends just on rapidity.
//
// Note that this has to be called BEFORE any attempt to do an
// actual computation
void GridMedianBackgroundEstimator::set_rescaling_class(const FunctionOfPseudoJet<double> * rescaling_class_in) {
  // The rescaling is taken into account when particles are set. So
  // you need to call set_particles again if you set the rescaling
  // class. We thus warn if there are already some available
  // particles
  if (_cache_available)
    _warning_rescaling.warn("GridMedianBackgroundEstimator::set_rescaling_class(): trying to set the rescaling class when there are already particles that have been set is dangerous: the rescaling will not affect the already existing particles resulting in mis-estimation of rho. You need to call set_particles() again before proceeding with any background estimation.");
  
  BackgroundEstimatorBase::set_rescaling_class(rescaling_class_in);
}



FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh
