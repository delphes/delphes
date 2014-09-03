
//FJSTARTHEADER
// $Id: ClusterSequenceAreaBase.cc 3433 2014-07-23 08:17:03Z salam $
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




#include "fastjet/ClusterSequenceAreaBase.hh"
#include <algorithm>

FASTJET_BEGIN_NAMESPACE

using namespace std;


/// allow for warnings
LimitedWarning ClusterSequenceAreaBase::_warnings;
LimitedWarning ClusterSequenceAreaBase::_warnings_zero_area;
LimitedWarning ClusterSequenceAreaBase::_warnings_empty_area;

//----------------------------------------------------------------------
/// return the total area, within the selector's range, that is free
/// of jets.
/// 
/// Calculate this as (range area) - \sum_{i in range} A_i
///
/// for ClusterSequences with explicit ghosts, assume that there will
/// never be any empty area, i.e. it is always filled in by pure
/// ghosts jets. This holds for seq.rec. algorithms
double ClusterSequenceAreaBase::empty_area(const Selector & selector) const {

  if (has_explicit_ghosts()) {return 0.0;}
  else { return empty_area_from_jets(inclusive_jets(0.0), selector);}

}

//----------------------------------------------------------------------
/// return the total area, within range, that is free of jets.
/// 
/// Calculate this as (range area) - \sum_{i in range} A_i
///
double ClusterSequenceAreaBase::empty_area_from_jets(
                      const std::vector<PseudoJet> & all_jets,
                      const Selector & selector) const {
  _check_selector_good_for_median(selector);

  double empty = selector.area();
  for (unsigned i = 0; i < all_jets.size(); i++) {
    if (selector.pass(all_jets[i])) empty -= area(all_jets[i]);
  }
  return empty;
}

double ClusterSequenceAreaBase::median_pt_per_unit_area(const Selector & selector) const {
  return median_pt_per_unit_something(selector,false);
}

double ClusterSequenceAreaBase::median_pt_per_unit_area_4vector(const Selector & selector) const {
  return median_pt_per_unit_something(selector,true);
}


//----------------------------------------------------------------------
/// the median of (pt/area) for jets contained within range, counting
/// the empty area as if it were made up of a collection of empty
/// jets each of area (0.55 * pi R^2).
double ClusterSequenceAreaBase::median_pt_per_unit_something(
                const Selector & selector, bool use_area_4vector) const {

  double median, sigma, mean_area;
  get_median_rho_and_sigma(selector, use_area_4vector, median, sigma, mean_area);
  return median;

}


//----------------------------------------------------------------------
/// fits a form pt_per_unit_area(y) = a + b*y^2 for jets in range. 
/// exclude_above allows one to exclude large values of pt/area from fit. 
/// use_area_4vector = true uses the 4vector areas.
void ClusterSequenceAreaBase::parabolic_pt_per_unit_area(
       double & a, double & b, const Selector & selector, 
       double exclude_above, bool use_area_4vector) const {
  // sanity check on the selector: we require a finite area and that
  // it applies jet by jet (see BackgroundEstimator for more advanced
  // usage)
  _check_selector_good_for_median(selector);

  int n=0;
  int n_excluded = 0;
  double mean_f=0, mean_x2=0, mean_x4=0, mean_fx2=0; 

  vector<PseudoJet> incl_jets = inclusive_jets();

  for (unsigned i = 0; i < incl_jets.size(); i++) {
    if (selector.pass(incl_jets[i])) {
      double this_area;
      if ( use_area_4vector ) {
          this_area = area_4vector(incl_jets[i]).perp();     
      } else {
          this_area = area(incl_jets[i]);
      }
      double f = incl_jets[i].perp()/this_area;
      if (exclude_above <= 0.0 || f < exclude_above) {
	double x = incl_jets[i].rap(); double x2 = x*x;
	mean_f   += f;
	mean_x2  += x2;
	mean_x4  += x2*x2;
	mean_fx2 += f*x2;
	n++;
      } else {
	n_excluded++;
      }
    }
  }

  if (n <= 1) {
    // meaningful results require at least two jets inside the
    // area -- mind you if there are empty jets we should be in 
    // any case doing something special...
    a = 0.0;
    b = 0.0;
  } else {
    mean_f   /= n;
    mean_x2  /= n;
    mean_x4  /= n;
    mean_fx2 /= n;
    
    b = (mean_f*mean_x2 - mean_fx2)/(mean_x2*mean_x2 - mean_x4);
    a = mean_f - b*mean_x2;
  }
  //cerr << "n_excluded = "<< n_excluded << endl;
}



void ClusterSequenceAreaBase::get_median_rho_and_sigma(
            const Selector & selector, bool use_area_4vector,
            double & median, double & sigma, double & mean_area) const {

  vector<PseudoJet> incl_jets = inclusive_jets();
  get_median_rho_and_sigma(incl_jets, selector, use_area_4vector,
			   median, sigma, mean_area, true);
}


void ClusterSequenceAreaBase::get_median_rho_and_sigma(
            const vector<PseudoJet> & all_jets,
            const Selector & selector, bool use_area_4vector,
            double & median, double & sigma, double & mean_area,
	    bool all_are_incl) const {

  _check_jet_alg_good_for_median();

  // sanity check on the selector: we require a finite area and that
  // it applies jet by jet (see BackgroundEstimator for more advanced
  // usage)
  _check_selector_good_for_median(selector);

  vector<double> pt_over_areas;
  double total_area  = 0.0;
  double total_njets = 0;

  for (unsigned i = 0; i < all_jets.size(); i++) {
    if (selector.pass(all_jets[i])) {
      double this_area;
      if (use_area_4vector) {
          this_area = area_4vector(all_jets[i]).perp();
      } else {
          this_area = area(all_jets[i]);
      }

      if (this_area>0) {
	pt_over_areas.push_back(all_jets[i].perp()/this_area);
      } else {
	_warnings_zero_area.warn("ClusterSequenceAreaBase::get_median_rho_and_sigma(...): discarded jet with zero area. Zero-area jets may be due to (i) too large a ghost area (ii) a jet being outside the ghost range (iii) the computation not being done using an appropriate algorithm (kt;C/A).");
      }

      total_area  += this_area;
      total_njets += 1.0;
    }
  }

  // there is nothing inside our region, so answer will always be zero
  if (pt_over_areas.size() == 0) {
    median = 0.0;
    sigma  = 0.0;
    mean_area = 0.0;
    return;
  }
  
  // get median (pt/area) [this is the "old" median definition. It considers
  // only the "real" jets in calculating the median, i.e. excluding the
  // only-ghost ones; it will be supplemented with more info below]
  sort(pt_over_areas.begin(), pt_over_areas.end());

  // now get the median & error, accounting for empty jets
  // define the fractions of distribution at median, median-1sigma
  double posn[2] = {0.5, (1.0-0.6827)/2.0};
  double res[2];
  
  double n_empty, empty_a;
  if (has_explicit_ghosts()) {
    // NB: the following lines of code are potentially incorrect in cases
    //     where there are unclustered particles (empty_area would do a better job,
    //     at least for active areas). This is not an issue with kt or C/A, or other
    //     algorithms that cluster all particles (and the median estimation should in 
    //     any case only be done with kt or C/A!)
    empty_a = 0.0;
    n_empty = 0;
  } else if (all_are_incl) {
    // the default case
    empty_a = empty_area(selector);
    n_empty = n_empty_jets(selector);
  } else {
    // this one is intended to be used when e.g. one runs C/A, then looks at its
    // exclusive jets in order to get an effective smaller R value, and passes those
    // to this routine.
    empty_a = empty_area_from_jets(all_jets, selector);
    mean_area = total_area / total_njets; // temporary value
    n_empty   = empty_a / mean_area;
  }
  //cout << "*** tot_area = " << total_area << ", empty_a = " << empty_a << endl;
  //cout << "*** n_empty = " << n_empty << ", ntotal =  " << total_njets << endl;
  total_njets += n_empty;
  total_area  += empty_a;

  // we need an int (rather than an unsigned int) with the size of the
  // pt_over_areas array, because we'll often be doing subtraction of
  // -1, negating it, etc. All of these operations go crazy with unsigned ints.
  int pt_over_areas_size = pt_over_areas.size();
  if (n_empty < -pt_over_areas_size/4.0)
    _warnings_empty_area.warn("ClusterSequenceAreaBase::get_median_rho_and_sigma(...): the estimated empty area is suspiciously large and negative and may lead to an over-estimation of rho. This may be due to (i) a rare statistical fluctuation or (ii) too small a range used to estimate the background properties.");

  for (int i = 0; i < 2; i++) {
    double nj_median_pos = 
      (pt_over_areas_size-1.0 + n_empty)*posn[i] - n_empty;
    double nj_median_ratio;
    if (nj_median_pos >= 0 && pt_over_areas_size > 1) {
      int int_nj_median = int(nj_median_pos);
 
     // avoid potential overflow issues
      if (int_nj_median+1 > pt_over_areas_size-1){
	int_nj_median = pt_over_areas_size-2;
	nj_median_pos = pt_over_areas_size-1;
      }

      nj_median_ratio = 
        pt_over_areas[int_nj_median] * (int_nj_median+1-nj_median_pos)
        + pt_over_areas[int_nj_median+1] * (nj_median_pos - int_nj_median);
    } else {
      nj_median_ratio = 0.0;
    }
    res[i] = nj_median_ratio;
  }
  median = res[0];
  double error  = res[0] - res[1];
  mean_area = total_area / total_njets;
  sigma  = error * sqrt(mean_area);
}


/// return a vector of all subtracted jets, using area_4vector, given rho.
/// Only inclusive_jets above ptmin are subtracted and returned.
/// the ordering is the same as that of sorted_by_pt(cs.inclusive_jets()),
/// i.e. not necessarily ordered in pt once subtracted
vector<PseudoJet> ClusterSequenceAreaBase::subtracted_jets(const double rho,
                                                           const double ptmin) 
                                                           const {
  vector<PseudoJet> sub_jets;
  vector<PseudoJet> jets_local = sorted_by_pt(inclusive_jets(ptmin));
  for (unsigned i=0; i<jets_local.size(); i++) {
     PseudoJet sub_jet = subtracted_jet(jets_local[i],rho);
     sub_jets.push_back(sub_jet);
  }
  return sub_jets;
}

/// return a vector of subtracted jets, using area_4vector.
/// Only inclusive_jets above ptmin are subtracted and returned.
/// the ordering is the same as that of sorted_by_pt(cs.inclusive_jets()),
/// i.e. not necessarily ordered in pt once subtracted
vector<PseudoJet> ClusterSequenceAreaBase::subtracted_jets(
                                                 const Selector & selector, 
						 const double ptmin)
						 const {
  double rho = median_pt_per_unit_area_4vector(selector);
  return subtracted_jets(rho,ptmin);
}


/// return a subtracted jet, using area_4vector, given rho
PseudoJet ClusterSequenceAreaBase::subtracted_jet(const PseudoJet & jet,
                                                  const double rho) const {
  PseudoJet area4vect = area_4vector(jet);
  PseudoJet sub_jet;
  // sanity check
  if (rho*area4vect.perp() < jet.perp() ) { 
    sub_jet = jet - rho*area4vect;
  } else { sub_jet = PseudoJet(0.0,0.0,0.0,0.0); }
  
  // make sure the subtracted jet has the same index (cluster, user, csw)
  // (i.e. "looks like") the original jet
  sub_jet.set_cluster_hist_index(jet.cluster_hist_index());
  sub_jet.set_user_index(jet.user_index());
  // do not use CS::_set_structure_shared_ptr here, which should
  // only be called to maintain the tally during construction
  sub_jet.set_structure_shared_ptr(jet.structure_shared_ptr());
  return sub_jet;
}


/// return a subtracted jet, using area_4vector;  note that this is
/// potentially inefficient if repeatedly used for many different
/// jets, because rho will be recalculated each time around.
PseudoJet ClusterSequenceAreaBase::subtracted_jet(const PseudoJet & jet,
                                       const Selector & selector) const {
  double rho = median_pt_per_unit_area_4vector(selector);
  PseudoJet sub_jet = subtracted_jet(jet, rho);
  return sub_jet;
}


/// return the subtracted pt, given rho
double ClusterSequenceAreaBase::subtracted_pt(const PseudoJet & jet,
                                              const double rho,
                                              bool use_area_4vector) const {
  if ( use_area_4vector ) { 
     PseudoJet sub_jet = subtracted_jet(jet,rho);
     return sub_jet.perp();
  } else {
     return jet.perp() - rho*area(jet);
  }
}  


/// return the subtracted pt; note that this is
/// potentially inefficient if repeatedly used for many different
/// jets, because rho will be recalculated each time around.
double ClusterSequenceAreaBase::subtracted_pt(const PseudoJet & jet,
                                              const Selector & selector,
                                              bool use_area_4vector) const {
  if ( use_area_4vector ) { 
     PseudoJet sub_jet = subtracted_jet(jet,selector);
     return sub_jet.perp();
  } else {
     double rho = median_pt_per_unit_area(selector);
     return subtracted_pt(jet,rho,false);
  }
}  

// check the selector is suited for the computations i.e. applies jet
// by jet and has a finite area
void ClusterSequenceAreaBase::_check_selector_good_for_median(const Selector &selector) const{
  // make sure the selector has a finite area
  if ((! has_explicit_ghosts()) &&  (! selector.has_finite_area())){
    throw Error("ClusterSequenceAreaBase: empty area can only be computed from selectors with a finite area");
  }

  // make sure the selector applies jet by jet
  if (! selector.applies_jet_by_jet()){
    throw Error("ClusterSequenceAreaBase: empty area can only be computed from selectors that apply jet by jet");
  }
}


/// check the jet algorithm is suitable (and if not issue a warning)
void ClusterSequenceAreaBase::_check_jet_alg_good_for_median() const {
  if (jet_def().jet_algorithm() != kt_algorithm
      && jet_def().jet_algorithm() != cambridge_algorithm
      && jet_def().jet_algorithm() !=  cambridge_for_passive_algorithm) {
    _warnings.warn("ClusterSequenceAreaBase: jet_def being used may not be suitable for estimating diffuse backgrounds (good options are kt, cam)");
  }
}



FASTJET_END_NAMESPACE
