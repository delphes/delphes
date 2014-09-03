//FJSTARTHEADER
// $Id: BackgroundEstimatorBase.cc 3433 2014-07-23 08:17:03Z salam $
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


#include "fastjet/tools/BackgroundEstimatorBase.hh"

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

LimitedWarning BackgroundEstimatorBase::_warnings_empty_area;

//----------------------------------------------------------------------
// given a quantity in a vector (e.g. pt_over_area) and knowledge
// about the number of empty jets, calculate the median and
// stand_dev_if_gaussian (roughly from the 16th percentile)
//
// If do_fj2_calculation is set to true then this performs FastJet
// 2.X estimation of the standard deviation, which has a spurious
// offset in the limit of a small number of jets.
void BackgroundEstimatorBase::_median_and_stddev(const vector<double> & quantity_vector, 
                                                 double n_empty_jets, 
                                                 double & median, 
                                                 double & stand_dev_if_gaussian,
                                                 bool do_fj2_calculation) const {

  // this check is redundant (the code below behaves sensibly even
  // with a zero size), but serves as a reminder of what happens if
  // the quantity vector is zero-sized
  if (quantity_vector.size() == 0) {
    median = 0;
    stand_dev_if_gaussian = 0;
    return;
  }

  vector<double> sorted_quantity_vector = quantity_vector;
  sort(sorted_quantity_vector.begin(), sorted_quantity_vector.end());

  // empty area can sometimes be negative; with small ranges this can
  // become pathological, so warn the user
  int n_jets_used = sorted_quantity_vector.size();
  if (n_empty_jets < -n_jets_used/4.0)
    _warnings_empty_area.warn("BackgroundEstimatorBase::_median_and_stddev(...): the estimated empty area is suspiciously large and negative and may lead to an over-estimation of rho. This may be due to (i) a rare statistical fluctuation or (ii) too small a range used to estimate the background properties.");

  // now get the median & error, accounting for empty jets;
  // define the fractions of distribution at median, median-1sigma
  double posn[2] = {0.5, (1.0-0.6827)/2.0};
  double res[2];
  for (int i = 0; i < 2; i++) {
    res[i] = _percentile(sorted_quantity_vector, posn[i], n_empty_jets, 
                         do_fj2_calculation);
  }
  
  median = res[0];
  stand_dev_if_gaussian = res[0] - res[1];
}


//----------------------------------------------------------------------
// computes a percentile of a given _sorted_ vector of quantities
//  - sorted_quantities        the (sorted) vector contains the data sample
//  - percentile               the percentile (defined between 0 and 1) to compute
//  - nempty                   an additional number of 0's
//                             (considered at the beginning of 
//                             the quantity vector)
//  - do_fj2_calculation       carry out the calculation as it
//                             was done in fj2 (suffers from "edge effects")
double BackgroundEstimatorBase::_percentile(const vector<double> & sorted_quantities, 
                                            const double percentile, 
                                            const double nempty,
                                            const bool do_fj2_calculation
                                            ) const {
  assert(percentile >= 0.0 && percentile <= 1.0);

  int quantities_size = sorted_quantities.size();
  if (quantities_size == 0) return 0;

  double total_njets = quantities_size + nempty;
  double percentile_pos;
  if (do_fj2_calculation) {
    percentile_pos = (total_njets-1)*percentile - nempty;
  } else {
    percentile_pos = (total_njets)*percentile - nempty - 0.5;
  }

  double result;
  if (percentile_pos >= 0 && quantities_size > 1) {
    int int_percentile_pos = int(percentile_pos);

    // avoid potential overflow issues
    if (int_percentile_pos+1 > quantities_size-1){
      int_percentile_pos = quantities_size-2;
      percentile_pos = quantities_size-1;
    }

    result =
      sorted_quantities[int_percentile_pos] * (int_percentile_pos+1-percentile_pos)
      + sorted_quantities[int_percentile_pos+1] * (percentile_pos - int_percentile_pos);
    

  } else if (percentile_pos > -0.5 && quantities_size >= 1 
             && !do_fj2_calculation) {
    // in the LHS of this "bin", just keep a constant value (we could have
    // interpolated to zero, but this might misbehave in cases where all jets
    // are active, because it would go to zero too fast)
    result = sorted_quantities[0];
  } else {
    result = 0.0;
  }
  return result;


}


FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh
